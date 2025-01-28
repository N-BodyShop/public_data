# N-Body Shop Public Data Release Tools

The paper is started in [Overleaf](https://www.overleaf.com/project/66969b53942ad127acd619db) see @trquinn if you want to edit.

## TANGOS Database
The simplest way to build a standard, N-Body Shop TANGOS database is to run the
`build_tangos_DB.sh PATH_TO_SIMULATIONS PATH_TO_DB` script.  Your
`PATH_TO_SIMULATIONS` directory should contain one directory for each simulation
you want packaged in the database, with these directories containing your
tipsy/N-Chilada outputs.  For example, if I wanted to add all of the MUGS2
simulations stored in `/home/kellerbw/data/MUGS2` to a TANGOS database stored at
`/home/kellerbw/MUGS2.db`, I would use:

```
./build_tangos_DB.sh /home/kellerbw/data/MUGS2 /home/kellerbw/MUGS2.db
```

This script _MUST_ be passed absolute paths: no wildcards or tilde-prefixes.
