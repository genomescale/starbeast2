# For developer

These scripts are used in the developer env, not for user.

## Usage

```
bin/beast -validate path/to/file.xml     # structural check only, no MCMC run
bin/beast -overwrite path/to/file.xml    # full MCMC run, overwrite existing logs
```

`bin/beast` delegates to the `exec-maven-plugin` already configured in
`pom.xml` (`mvn compile exec:exec -Dbeast.args=...`), so Maven resolves the
full module path (beast-base, beast-fx, SA, MM, guava, commons-math3, ...)
straight from `pom.xml`'s own dependency graph. 
