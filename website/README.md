# DELFIN website draft

Local review only. Nothing has been published or configured on GitHub Pages.

## Preview

Open `index.html` in a browser directly, or serve this folder from the repository:

```sh
python -m http.server 8080 --bind 127.0.0.1 --directory website
```

Then visit `http://localhost:8080` **on the same computer**. If the repository is
on a remote machine, forward the port with SSH, or copy the `website` folder
to your computer and open `index.html`. No installation or build step is needed.

## Change and extend

- `content.js`: hero text, workflow stages, research areas, Animal Loop and videos.
- `index.html`: section order, mission text, interface descriptions, installation,
  metadata and links. New sections can be added without changing the others.
- `styles.css`: design tokens, responsive layouts and motion. Brand colours are
  at the beginning. Existing logos are copied unchanged into `assets/`.
- `app.js`: canvas illustration, selectors, scroll reveals, motion control,
  copy button and video dialog. No external runtime dependencies.
- `PLAN.md`: concept and source material behind the draft.

All page copy is editable draft copy. The site is English in this first version.
The molecular viewer displays actual MANTA-generated fac-Ir(ppy)3 coordinates;
see `assets/STRUCTURE.md` for the reproducible command and provenance. It is a
starting geometry, not a DFT result. The terminal sequence remains illustrative.
WEDDELL is included as part of the Animal Loop without an invented capability.
Personal profile teasers and links have been removed at the user's request.

The second visual direction uses a dark scientific structure explorer, sans-serif
editorial headings and restrained motion. The original blue identity, mission and
Animal Loop remain. The DELFIN wordmark has no trailing full stop.

### Add videos

Add an entry to `videos` in `content.js`. The section appears automatically when
the array is nonempty. A local video example:

```js
videos: [
  {
    title: 'Your first DELFIN calculation',
    description: 'From molecular input to a documented result.',
    type: 'file',
    src: 'assets/videos/first-calculation.mp4',
    poster: 'assets/videos/first-calculation.webp',
    captions: 'assets/videos/first-calculation.en.vtt',
    language: 'en',
    captionLabel: 'English'
  }
]
```

Place the actual files at those paths. Use a local poster and captions. For large
videos, use an external video host rather than storing recordings in Git history.
YouTube entries use `type: 'youtube'` and
`src: 'https://www.youtube-nocookie.com/embed/VIDEO_ID'`, with a local poster.
The player is created only after a click; closing it stops playback and removes
the iframe. YouTube captions are managed in YouTube. The site itself does not
autoplay videos or contact YouTube before the visitor opens one.

### Later growth

Add tutorial pages, a publications list or a changelog as ordinary HTML pages,
using the shared stylesheet. No need to replace the landing page. Keep media
paths relative so the site works at `/DELFIN/` on GitHub Pages. If the volume of
documentation grows substantially, a separate documentation generator can be
introduced later without changing the landing page’s visual identity.

## Review and release

Review text, applications, the Animal Loop, motion and mobile presentation with
Max before release. After approval, a Pages workflow can publish **only this
folder**. No deployment workflow is enabled in this draft.
