import chembl_downloader


def main():
    path = chembl_downloader.download_extract_sqlite(version="35")
    print(path)


if __name__ == "__main__":
    main()
