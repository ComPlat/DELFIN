# Ratings

`ratings/<set>/review_<name>.json` -- the ratings of one reviewer for the
review set `review_sets/<set>/`, as written by the Review tab and sent with
`delfin review submit` (one file per reviewer; a pull request adds or
replaces only its own file). `delfin review merge` brings them together and
measures the agreement; see [docs/REVIEW_TEAM_GUIDE.md](../docs/REVIEW_TEAM_GUIDE.md).
