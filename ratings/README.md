# Ratings

`ratings/<set>/review_<name>.json` -- the ratings of one reviewer for the
review set `review_sets/<set>/`, as written by the Review tab and sent with
`delfin review submit` (one file per reviewer; a pull request adds or
replaces only its own file). `ratings/<set>/claims/<name>.json` -- a
reviewer's reserved block (`delfin review claim`; it runs out after 28 days).
`ratings/<set>/owners/<name>.json` -- the GitHub account a name belongs to,
added by the first pull request under that name.
`delfin review status` and `delfin review merge` bring them together and
measure the agreement; every pull request into this folder is checked by
`.github/workflows/ratings.yml`. See
[docs/REVIEW_TEAM_GUIDE.md](../docs/REVIEW_TEAM_GUIDE.md).
