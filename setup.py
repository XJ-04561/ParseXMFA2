from setuptools import setup, find_packages

setup(
    name="ParseXMFA2",
    version=1.0,
    url=None,
    description="XMFA parser using `mmap` to index the blocks in order to look up SNPs and export results to a .vcf format..",

    # Author details
    author="Fredrik Sörensen",
    author_email="fredrik.sorensen@foi.se",

    license='GNU GENERAL PUBLIC LICENSE version 3',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.8",

    keywords="variant-calling snp-calling",

    install_requires=["VariantCallFixer"],
    entry_points={"console_scripts": [
                    "ParseXMFA2=ParseXMFA2.ParseXMFA2"
    ]})
