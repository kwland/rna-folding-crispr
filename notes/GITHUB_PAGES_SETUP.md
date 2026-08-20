# GitHub Pages Setup

The website is in the `docs/` folder so it can be published directly with GitHub Pages.

## Publish

1. Push this repository to GitHub.
2. Open the repository on GitHub.
3. Go to **Settings** -> **Pages**.
4. Under **Build and deployment**, choose **Deploy from a branch**.
5. Select your main branch and set the folder to `/docs`.
6. Save. GitHub will publish the site after the Pages build finishes.

## Local Preview

From the project folder:

```bash
python -m http.server 8000 -d docs
```

Then open:

```text
http://127.0.0.1:8000/
```

The site is static: it does not need Flask, Node, React, or a database.
