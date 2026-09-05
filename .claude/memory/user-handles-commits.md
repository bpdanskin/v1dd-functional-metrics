---
name: user-handles-commits
description: "Leave changes in the working tree — the user inspects the diff and commits themselves."
metadata:
  node_type: memory
  type: feedback
  modified: 2026-09-03
---

**Do not run `git commit`.** Leave finished work as uncommitted changes and say what
changed and in which files.

**Why:** the user inspects diffs before committing, and that review is part of how they
stay oriented in a codebase an agent is changing quickly. Committing for them skips it.

**How to apply:** finish the work, summarise what changed and where, and stop. If a commit
is genuinely needed, ask. A one-off request — "help me with this merge" — is permission for
*that* commit only, not a standing grant.
