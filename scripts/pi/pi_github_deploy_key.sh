#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Give the agent push access to exactly one GitHub repository via a deploy key.
#
# Runs ON the Pi as the agent user. Idempotent. Alternative to a personal
# token when the organisation does not allow fine-grained tokens: a deploy
# key reaches one repository only and is added by any repository admin.
#
#   1. create ~/.ssh/delfin_github_deploy (ed25519) unless it exists
#   2. pin github.com's host key, only if its fingerprint is GitHub's
#      published one
#   3. use the key for github.com and switch the checkout's origin to SSH
#   4. print the public key and where to add it ("Allow write access")
#
# Usage:  pi_github_deploy_key.sh [owner/repo] [checkout]
#         defaults: ComPlat/DELFIN  ~/delfin/work/DELFIN
# Then:   pi_github_deploy_key.sh --check   (auth + push permission, dry run)
# ============================================================================

# https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/githubs-ssh-key-fingerprints
GITHUB_ED25519_FP="SHA256:+DiY3wvvV6TuJJhbpZisF/zLDA0zPMSvHdkr4UvCOqU"
KEY="$HOME/.ssh/delfin_github_deploy"

check=0
if [[ "${1:-}" == "--check" ]]; then check=1; shift; fi
REPO="${1:-ComPlat/DELFIN}"
CHECKOUT="${2:-$HOME/delfin/work/DELFIN}"

if [[ $check == 1 ]]; then
    ssh -n -T git@github.com 2>&1 | head -1 || true
    git -C "$CHECKOUT" push --dry-run origin "HEAD:refs/heads/delfin-write-probe" </dev/null
    exit 0
fi

install -d -m 700 "$HOME/.ssh"
[[ -f "$KEY" ]] || ssh-keygen -q -t ed25519 -N "" -C "$(hostname) deploy key ($REPO)" -f "$KEY"

scan="$(ssh-keyscan -t ed25519 github.com 2>/dev/null)"
fp="$(ssh-keygen -lf - <<<"$scan" | awk '{print $2}')"
if [[ "$fp" != "$GITHUB_ED25519_FP" ]]; then
    echo "error: github.com host key $fp is not GitHub's published $GITHUB_ED25519_FP" >&2
    exit 1
fi
touch "$HOME/.ssh/known_hosts"
grep -qF "$(cut -d' ' -f2- <<<"$scan")" "$HOME/.ssh/known_hosts" || echo "$scan" >> "$HOME/.ssh/known_hosts"

if ! grep -q "^Host github.com" "$HOME/.ssh/config" 2>/dev/null; then
    cat >> "$HOME/.ssh/config" <<EOF
Host github.com
    User git
    IdentityFile $KEY
    IdentitiesOnly yes
EOF
fi
chmod 600 "$HOME/.ssh/config"

if [[ -d "$CHECKOUT/.git" ]]; then
    git -C "$CHECKOUT" remote set-url origin "git@github.com:$REPO.git"
fi

echo "Add this key at https://github.com/$REPO/settings/keys/new"
echo "with \"Allow write access\" ticked:"
echo
cat "$KEY.pub"
echo
echo "Then: $0 --check"
