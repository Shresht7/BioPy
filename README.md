# 🧬 Learning BioInformatics

Python and Bioinformatics fundamentals.

---

## 📕 References

- [rosalind.info](https://rosalind.info/problems/locations/)
- [YouTube Playlist: BioInformatics in Python - rebelScience](https://www.youtube.com/playlist?list=PLpSOMAcxEB_jUKMvdl8rHqNiZXFIrtd5G)
- [YouTube: Bioinformatics with Biopython - Lana Dominkovic](https://www.youtube.com/watch?v=ocA2IMe7dpA)

---

## Code

### Virtual Environment

#### To create a virtual environment

```sh
python -m venv .venv
```

This will create a `.venv` folder. If it already exists, you can move onto the next step.

#### Active the virtual environment

Run the appropriate shell script in the `.venv` folder

```sh
. ./.venv/bin/activate
```

or 

```pwsh
. .\.venv\Scripts\Activate.ps1
```

#### Install dependencies

To install all requirements

```sh
pip install --upgrade -r requirements.txt
```

>[!NOTE]
> To save the dependencies run: `pip freeze > requirements.txt`

---

## License

This project is licensed under the MIT License. See [LICENSE](./LICENSE) file for more details.
