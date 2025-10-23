# 🧬 Learning BioInformatics

Python and Bioinformatics fundamentals.

> [!CAUTION]
> **DO NOT** take any of this work seriously. I am inept at Biology, and asked LLMs for help. Proceed with extreme skepticism.

---

## 📕 References

- [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/)
- [RCSB Protein Data Bank](https://www.rcsb.org/)
- [BioPython Documentation](https://biopython.org/)
- [rosalind.info](https://rosalind.info/problems/locations/)
- [YouTube: The Structure of DNA - MITx Bio](https://www.youtube.com/watch?v=o_-6JXLYS-k)
- [YouTube: DNA Structure and Replication: Crash Course Biology #10](https://www.youtube.com/watch?v=8kK2zwjRV0M)
- [YouTube: DNA Transcription & Translation: Crash Course Biology #11](https://www.youtube.com/watch?v=itsb2SqR-R0&t=52s)
- [YouTube Playlist: BioInformatics in Python - rebelScience](https://www.youtube.com/playlist?list=PLpSOMAcxEB_jUKMvdl8rHqNiZXFIrtd5G)
- [YouTube: Bioinformatics with Biopython - Lana Dominkovic](https://www.youtube.com/watch?v=ocA2IMe7dpA)
- [GitHub: 12 Days of BioPython - Lana Dominkovic](https://github.com/lanadominkovic/12-days-of-biopython)
- [YouTube: Introduction to Computational and Systems Biology - MIT OpenCourseWare](https://www.youtube.com/watch?v=lJzybEXmIj0)

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
