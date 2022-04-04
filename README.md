# ftol_data

Home of [Fern Tree of Life (FTOL)](https://fernphy.github.io) phylogeny and associated data files.

The data in this repo are automatically generated by the `ftol` repo (https://github.com/fernphy/ftol).

These files should not be edited by hand. Rather, the [`snapshot_ftol_data.R`](https://github.com/fernphy/ftol/blob/main/R/snapshot_ftol_data.R) script should be run to commit changes and automatically push to the [remote repository](https://github.com/fernphy/ftol_data).

For more information, please see the accompanying paper:
- Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. An open and continuously updated fern tree of life (FTOL) https://doi.org/10.1101/2022.03.31.486640 

## Commit messages

The commit messages in this repo have a special format. They have two parts, `code` and `comment`.

- The `code` part is formatted like `code=...`, where the part after the equals sign is a commit hash in the [ftol](https://github.com/fernphy/ftol) repo.
- The `comment` part is formatted like `comment=...`, where the part after the equals sign is coment for that commit in the [ftol](https://github.com/fernphy/ftol) repo.

The files contained here are generated by code at the corresponding commit in the [ftol](https://github.com/fernphy/ftol) repo.

For example, [files with the commit](https://github.com/fernphy/ftol_data/commit/9051467d88bf606d897a799b5c1cce367fef4e42)

```
code=b54fffbb7dd5d4e959d339889b888f031202f7c6
comment=Update ftol_data_readme
```

were generated by [this code](https://github.com/fernphy/ftol/commit/b54fffbb7dd5d4e959d339889b888f031202f7c6).

## License

- [CC0](LICENSE)
