# How to create Sphinx docs

With the dev environment activated, and Sphinx installed you can create the html version by running the following command from this `docs` directory:

```bash
make html
```

And for the pdf version use

```bash
make latexpdf
```

Note: this last command requires a latex installation, which Jasmin servers don't seem to have.

```bash
open _build/html/index.html 
```

## Other important commands

To update the module references in the rst files

```bash
sphinx-apidoc -f -o . ..
```

```bash
sphinx-apidoc -f -o . ../ocean/SRC/
```

## Symbolic links

`docs/README.md` and `docs/gifs/`

are both symbolic links to items in the main directory.

This was done to trick sphinx into working, and seems to have worked so far.
