# Publishing this folder to GitHub

Recommended repository name: `clear-onyx-exact-kernel`

Suggested description:

> Five auditable mathematical papers and the Clear_Onyx computational kernel exposition, with LaTeX sources, checksums, regression witnesses, and continuous verification.

Suggested topics:

`mathematics`, `linear-algebra`, `dynamical-systems`, `projective-geometry`, `representation-theory`, `computational-mathematics`, `latex`, `reproducible-research`

## Browser route

1. Create a new empty repository on GitHub.
2. Do not initialize it with a README, `.gitignore`, or license; those files are already present here.
3. Open **Add file → Upload files**.
4. Drag the contents of this folder into the upload area.
5. Commit with: `Publish five-paper set and Clear_Onyx exposition`.
6. Open the **Actions** tab and confirm that **Verify release** passes.

## GitHub CLI route

From this folder:

```bash
git init -b main
git add .
git commit -m "Publish five-paper set and Clear_Onyx exposition"
gh repo create clear-onyx-exact-kernel --public --source=. --remote=origin --push
```

Before the commit, Git may ask you to set your author name and email. Use the identity you want displayed in the repository history.
