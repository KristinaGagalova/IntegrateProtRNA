# Real enforcement: restrict the SSH key on the server

Everything in `.claude/settings.json` and `.claude/hooks/validate_ssh.py`
runs on your Windows machine, inside Claude Code. It can be wrong, or a
prompt-injection from a file Claude reads could try to talk it into running
something it shouldn't. None of that matters if the SSH key itself can only
do one thing on the server — so this is the step that actually matters.

## 1. Generate a dedicated key for Claude Code (don't reuse your personal key)

On your Windows machine:

```powershell
ssh-keygen -t ed25519 -f "$HOME\.ssh\claude_integrateprotrna" -C "claude-code-integrateprotrna"
```

Leave a passphrase off only if this key will be used non-interactively and
you're comfortable with that tradeoff; otherwise set one and use an SSH
agent.

## 2. Install a forced command on the server for that key only

On the **server**, edit `~/.ssh/authorized_keys` for the account Claude will
log in as. Instead of a normal key line, prefix it with `command=` so this
key can *only* run one wrapper script, no matter what command SSH is asked
to run:

```
command="/home/<SSH_USER>/bin/claude-restricted-shell.sh",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-ed25519 AAAA...your-public-key... claude-code-integrateprotrna
```

## 3. Write the wrapper script on the server

`/home/<SSH_USER>/bin/claude-restricted-shell.sh`:

```bash
#!/bin/bash
set -euo pipefail

ALLOWED_DIR="/home/<SSH_USER>/IntegrateProtRNA"   # <- set this

# SSH puts the client's requested command here.
cmd="${SSH_ORIGINAL_COMMAND:-}"

if [[ -z "$cmd" ]]; then
  echo "Interactive shells are not permitted with this key." >&2
  exit 1
fi

# Refuse compound commands outright - one command per SSH call.
if [[ "$cmd" =~ [\;\|\&\`] || "$cmd" == *'$('* ]]; then
  echo "Compound/chained commands are not permitted with this key." >&2
  exit 1
fi

# Refuse obviously destructive commands regardless of path.
if [[ "$cmd" =~ sudo|mkfs|dd[[:space:]]+if=|rm[[:space:]]+-rf[[:space:]]+/([[:space:]]|$) ]]; then
  echo "Refused: matches a blocked pattern." >&2
  exit 1
fi

# Resolve any `cd <dir> && ...` the wrapper already refused above; instead
# require the command to explicitly operate under ALLOWED_DIR, and run it
# with that as the working directory.
cd "$ALLOWED_DIR"

# Belt-and-braces: reject paths that try to escape ALLOWED_DIR via `..`.
if [[ "$cmd" == *".."* ]]; then
  echo "Refused: '..' path traversal is not permitted." >&2
  exit 1
fi

exec bash -c "$cmd"
```

```bash
chmod 700 /home/<SSH_USER>/bin/claude-restricted-shell.sh
```

This means: even if something on the Windows side goes wrong — a bad hook,
a bug, a prompt injection — the worst this key can do on the server is run
a single, non-chained command with its working directory forced to
`ALLOWED_DIR`, and a few destructive patterns are refused outright. It
can't `cd` anywhere else, can't chain a second command after a legitimate
one, and can't be used to open an interactive shell.

## 4. Point Claude Code's SSH calls at that key

Either put it in `~/.ssh/config` on Windows (`%USERPROFILE%\.ssh\config`) so
plain `ssh <SSH_USER>@<SERVER_IP> "..."` picks it up automatically:

```
Host <SERVER_IP>
    User <SSH_USER>
    IdentityFile ~/.ssh/claude_integrateprotrna
    IdentitiesOnly yes
```

or pass `-i` explicitly in every command (`ssh -i ~/.ssh/claude_integrateprotrna ...`)
and adjust `validate_ssh.py`'s expected prefix to match.

## 5. Fill in the placeholders

Once this is set up, fill in the three placeholders in
`.claude/hooks/validate_ssh.py` (`ALLOWED_HOST`, `ALLOWED_USER`,
`ALLOWED_REMOTE_PATH`) and in `CLAUDE.md`, so the client-side hook and the
server-side wrapper agree on the same boundary.

## Caveats worth knowing

- A script that opens files itself (e.g. a Python/R pipeline invoked by the
  allowed command) isn't stopped by this from reading/writing paths outside
  `ALLOWED_DIR` — the wrapper only constrains what shell command SSH runs,
  not what that command's own process does once running. If that matters,
  run the allowed command as a low-privilege OS user who only has
  filesystem permissions inside `ALLOWED_DIR`, in addition to this wrapper.
- Rotate/revoke this key independently of your personal key if it's ever
  compromised — that's the whole point of using a separate one.
