/**
 * A credential's VALUE masked in a command, the shape of the command kept:
 * the REPL's rule (cli_approve.redact), so the window shows a command the
 * way chat shows it. `export TOKEN=sk-live-1234 && curl x` becomes
 * `export TOKEN=<redacted> && curl x`; a line that looks like a credential
 * with no assignment to mask is withheld whole. The worker masks commands
 * too, before they are sent or logged; this covers what reaches the window
 * from elsewhere.
 */

/** A key or filename that matches this is a credential until proven otherwise. */
const CRED_KEY = /(api[_-]?key|token|secret|password|passwd|credential|auth|\.env|\.pem|\.key|\.netrc|cookie|session[_-]?id)/i

/** NAME=value and --flag value after a credential-looking name. */
const CRED_ASSIGNMENT = /((?:api[_-]?key|token|secret|password|passwd|credential)[\w.-]*\s*=\s*|--(?:api[_-]?key|token|secret|password)[= ]\s*)(['"]?)([^\s'";|&]+)\2/gi

export function redact(value: string): string {
  if (!value) return value
  let hits = 0
  const masked = value.replace(CRED_ASSIGNMENT, (_match, lead: string) => {
    hits++
    return `${lead}<redacted>`
  })
  if (hits) return masked
  if (CRED_KEY.test(value)) return "<redacted: looks like a credential>"
  return value
}
