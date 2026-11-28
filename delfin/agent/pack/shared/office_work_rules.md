Nothing here fails loudly. A wrong number raises no exception, a letter
with a placeholder still in it gets sent, a value in the wrong field
looks perfectly correct. So the checking has to be deliberate:

1. **Read before you write.** Open the file and look at it. Never edit
   a cell reference you have not seen, never fill a field name you have
   not listed.
2. **Verify after you write.** Read the result back and compare it to
   what you intended. Say in your answer that you did, and what came
   back.
3. **Report coverage, not just results.** A total is a number plus the
   rows it covers: say which you read and name the ones skipped. A figure
   that silently dropped six rows is worse than no figure.
