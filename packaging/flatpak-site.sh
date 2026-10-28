#!/usr/bin/env bash
# Add the Flatpak repository to the Pages site (pages.yml): an OSTree repository at
# <site>/flatpak/repo holding the .flatpak bundle of the latest published release, plus
# fieldcraft.flatpakrepo, fieldcraft.flatpakref and an index page that says how to install.
#
# Usage: packaging/flatpak-site.sh SITE
#
# Environment:
#   GITHUB_REPOSITORY, GH_TOKEN  the repository whose releases are read (gh)
#   PAGES_BASE_URL               where SITE is served, e.g. https://tangency-labs.github.io/fieldcraft
#   RELEASE_TAG                  the published release to add (default: the latest published
#                                release that is not a pre-release)
#   FLATPAK_GPG_KEY              optional: an ASCII-armored GPG secret key without a passphrase.
#                                With it the commits and the summary are signed and the
#                                .flatpakrepo/.flatpakref carry its public key; without it the
#                                repository is unsigned (HTTPS only), and clients add it with
#                                GPG verification off.
#   GITHUB_OUTPUT                when set, gets `tag=<the release added>` (empty when none)
#
# The bundle is checked against the release's SHA256SUMS.txt, whose signature
# SHA256SUMS.txt.minisig must verify with the release key the app is built with (PUBLIC_KEY in
# crates/update/src/lib.rs); an unsigned release is refused.
# Every deploy builds the repository afresh from that bundle: it holds that one release (a
# bundle's commit has no parent, so there is no history to keep), and `flatpak update` moves
# installs to it. Before the first release with a .flatpak there is nothing to add and SITE is
# left as it is.
# Needs flatpak, ostree, gh, minisign and, for signing, gpg.
# After SolidCraft's packaging/flatpak-site.sh (MIT OR Apache-2.0; see ATTRIBUTION.md).
set -euo pipefail
# shellcheck source=env.sh
. "$(dirname "${BASH_SOURCE[0]}")/env.sh"

site="${1:?usage: packaging/flatpak-site.sh SITE}"
base_url="${PAGES_BASE_URL:?PAGES_BASE_URL is not set}"
base_url="${base_url%/}"
arch=x86_64
title=FieldCraft
comment="FieldCraft releases"
repo_url="$base_url/flatpak/repo/"
output() { if [ -n "${GITHUB_OUTPUT:-}" ]; then echo "$1" >>"$GITHUB_OUTPUT"; fi; }
# Leave the site without a repository.
skip() { warn "$1: the site has no Flatpak repository"; output tag=; exit 0; }

# The release: given, or the latest published one (GitHub's "latest" is never a draft or a
# pre-release).
tag="${RELEASE_TAG:-}"
if [ -z "$tag" ]; then
  tag="$(gh api "repos/$GITHUB_REPOSITORY/releases/latest" --jq .tag_name 2>/dev/null || true)"
  [ -n "$tag" ] || skip "no published release yet"
fi
info="$(gh api "repos/$GITHUB_REPOSITORY/releases/tags/$tag" --jq '"\(.draft) \(.prerelease)"')"
if [ "$info" != "false false" ]; then
  echo "error: $tag is not a published release (draft, pre-release: $info)" >&2
  exit 1
fi
version="${tag#v}"
asset="fieldcraft-$version-linux-$arch.flatpak"

work="$(mktemp -d)"
# The gpg-agent that signing starts goes with the key.
trap 'if [ -d "$work/gnupg" ]; then GNUPGHOME="$work/gnupg" gpgconf --kill gpg-agent; fi; rm -rf "$work"' EXIT
gh release download "$tag" --repo "$GITHUB_REPOSITORY" --dir "$work" --pattern SHA256SUMS.txt --pattern SHA256SUMS.txt.minisig --pattern "$asset"
if [ ! -s "$work/$asset" ]; then
  # Releases from before the Flatpak have no bundle. Only the default (latest) may be one.
  if [ -n "${RELEASE_TAG:-}" ]; then
    echo "error: release $tag has no $asset" >&2
    exit 1
  fi
  skip "release $tag has no $asset"
fi
# The checksums come from the release key, not just from the release page.
public="$(sed -n 's/^pub const PUBLIC_KEY: &str = "\(.*\)";$/\1/p' "$ROOT/crates/update/src/lib.rs")"
[ -n "$public" ] || { echo "error: no PUBLIC_KEY in crates/update/src/lib.rs" >&2; exit 1; }
[ -s "$work/SHA256SUMS.txt.minisig" ] || { echo "error: release $tag has no SHA256SUMS.txt.minisig" >&2; exit 1; }
minisign -Vqm "$work/SHA256SUMS.txt" -x "$work/SHA256SUMS.txt.minisig" -P "$public"
(cd "$work" && awk -v f="$asset" '$2 == f || $2 == "*" f { print; n++ } END { exit n != 1 }' SHA256SUMS.txt | sha256sum --check -)

sign=()
if [ -n "${FLATPAK_GPG_KEY:-}" ]; then
  export GNUPGHOME="$work/gnupg"
  (umask 077 && mkdir -p "$GNUPGHOME")
  printf '%s\n' "$FLATPAK_GPG_KEY" | gpg --batch --quiet --import
  key="$(gpg --batch --list-secret-keys --with-colons | awk -F: '$1 == "fpr" { print $10; exit }')"
  [ -n "$key" ] || { echo "error: FLATPAK_GPG_KEY holds no secret key" >&2; exit 1; }
  sign=(--gpg-sign="$key" --gpg-homedir="$GNUPGHOME")
  gpg_key_line="GPGKey=$(gpg --batch --export "$key" | base64 -w0)"
  echo "signing with $key"
else
  warn "FLATPAK_GPG_KEY is not set: the Flatpak repository is not signed (docs/releasing.md)"
  gpg_key_line=""
fi

repo="$site/flatpak/repo"
rm -rf "$site/flatpak"
mkdir -p "$repo"
ostree --repo="$repo" init --mode=archive-z2
ref="app/$APP_ID/$arch/stable"
flatpak build-import-bundle "${sign[@]}" --no-update-summary "$repo" "$work/$asset"
ostree --repo="$repo" rev-parse "$ref" >/dev/null || { echo "error: $asset holds no $ref" >&2; exit 1; }
flatpak build-update-repo "${sign[@]}" --title="$title" --comment="$comment" \
  --homepage="$base_url/" --default-branch=stable "$repo"
test -s "$repo/summary"
echo "$asset: $ref $(ostree --repo="$repo" rev-parse "$ref")"

{
  echo "[Flatpak Repo]"
  echo "Title=$title"
  echo "Comment=$comment"
  echo "Url=$repo_url"
  echo "Homepage=$base_url/"
  if [ -n "$gpg_key_line" ]; then echo "$gpg_key_line"; fi
} >"$site/flatpak/fieldcraft.flatpakrepo"
{
  echo "[Flatpak Ref]"
  echo "Title=$title"
  echo "Name=$APP_ID"
  echo "Branch=stable"
  echo "Url=$repo_url"
  echo "SuggestRemoteName=fieldcraft"
  echo "Homepage=$base_url/"
  echo "IsRuntime=false"
  echo "RuntimeRepo=$FLATHUB"
  if [ -n "$gpg_key_line" ]; then echo "$gpg_key_line"; fi
} >"$site/flatpak/fieldcraft.flatpakref"
cat >"$site/flatpak/index.html" <<HTML
<!doctype html>
<html lang="en">
<head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>FieldCraft Flatpak repository</title></head>
<body>
<h1>FieldCraft Flatpak repository</h1>
<p>Latest release: <a href="https://github.com/$GITHUB_REPOSITORY/releases/tag/$tag">$tag</a> (Linux x86_64).</p>
<p>Install it, with the freedesktop runtime from Flathub, and keep it up to date:</p>
<pre>flatpak install --user $base_url/flatpak/fieldcraft.flatpakref
flatpak run $APP_ID
flatpak update --user</pre>
<p>Or add the repository first (<a href="fieldcraft.flatpakrepo">fieldcraft.flatpakrepo</a>):</p>
<pre>flatpak remote-add --user --if-not-exists flathub $FLATHUB
flatpak remote-add --user --if-not-exists fieldcraft $base_url/flatpak/fieldcraft.flatpakrepo
flatpak install --user fieldcraft $APP_ID</pre>
</body>
</html>
HTML
output "tag=$tag"
