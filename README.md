# Ben Chaddha — Academic Website (Hugo + PaperMod)

A clean, static Hugo site for research writing, technical notes, and project portfolio entries.

## Stack

- [Hugo](https://gohugo.io/) (extended)
- [PaperMod](https://github.com/adityatelange/hugo-PaperMod) theme (git submodule)
- GitHub Actions for build + deploy to GitHub Pages

## Local development

```bash
git clone <your-repo-url>
cd benchaddha.github.io
git submodule update --init --recursive
hugo server -D
```

Open <http://localhost:1313>.

## Site structure

- `content/_index.md` — home
- `content/about/` — bio + interests
- `content/writing/` — blog posts
- `content/notes/` — short technical notes (math supported)
- `content/projects/` — portfolio entries
- `content/cv/` — CV page (+ PDF placeholder path)
- `content/contact/` — contact details

## Creating content

### New writing post

```bash
hugo new content writing/my-new-post.md
```

### New technical note

```bash
hugo new content notes/my-note.md
```

Set `math: true` in front matter when you want KaTeX rendering.

### New project entry

```bash
hugo new content projects/my-project.md
```

## Configuration notes

- Edit `hugo.toml` to update `baseURL`, title, menus, and social links.
- `baseURL` is set for project pages by default (`https://<username>.github.io/ben-site/`).
  - For user/organization pages, set it to `https://<username>.github.io/`.
- RSS is enabled for home, sections, tags, and categories.
- Syntax highlighting is enabled via Hugo Chroma.

## Deployment (GitHub Pages)

The workflow at `.github/workflows/hugo-pages.yml`:

1. Checks out the repo with submodules
2. Installs Hugo extended
3. Builds the site using a GitHub Pages-compatible base URL
4. Uploads and deploys `public/` to GitHub Pages

After pushing to `main` or `master`, enable Pages in your repository settings (Build and deployment → GitHub Actions).
