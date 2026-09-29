# Getting a DOI for this repository — step by step

Zenodo archives a snapshot of a GitHub release and mints a DOI for it. Total time: about ten minutes.

## 1. Push this repository to GitHub

```
cd Tile_18O_Isotopes
git init
git add .
git commit -m "Data and code for Arfania et al. 2026, JGR Biogeosciences"
git branch -M main
git remote add origin https://github.com/Arfania-coder/Tile_18O_Isotopes.git
git push -u origin main
```

If the repository already exists on GitHub with old content, push into it with `git push -u origin main --force`
(this replaces the old contents). Make sure the repository is set to **Public** in Settings.

Check from a logged-out browser that https://github.com/Arfania-coder/Tile_18O_Isotopes resolves.

## 2. Link GitHub to Zenodo

1. Go to https://zenodo.org and log in with your GitHub account (Log in → GitHub).
2. Click your name (top right) → **GitHub**.
3. Find `Arfania-coder/Tile_18O_Isotopes` in the list and switch it **On**.
   (If it is not listed, click "Sync now".)

## 3. Create a release on GitHub

1. On the repository page, click **Releases** (right-hand column) → **Create a new release**.
2. Tag: `v1.0.0`. Title: `v1.0.0 — JGR Biogeosciences submission`.
3. Description: `Data and code as used for manuscript 2026JG010070.`
4. Click **Publish release**.

Zenodo picks the release up automatically within a minute or two and mints a DOI. `.zenodo.json` in this
repository pre-fills the title, authors, description, keywords and licence, so nothing needs typing.

## 4. Get the DOI

On Zenodo → **Upload** (or your GitHub page there), the new record shows a DOI of the form
`10.5281/zenodo.NNNNNNN`. Zenodo issues two: a **version DOI** (this release) and a **concept DOI**
(all versions). Use the **version DOI** in the paper.

## 5. Put the DOI in the manuscript

Replace `[DOI]` in the Data Availability Statement, the Software Availability Statement, and the two
reference entries Arfania et al. (2026a) and Arfania (2026b):

```
https://doi.org/10.5281/zenodo.NNNNNNN
```

Also paste it into the `README.md` citation block, then commit and push (this does not change the DOI).

## If you later change the data or code

Create a new release (`v1.0.1`, `v1.1.0` …). Zenodo mints a new version DOI automatically. The concept DOI
always resolves to the latest version.
