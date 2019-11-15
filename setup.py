import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="i-impute",
    version="1.0.4",
    author="Xikang Feng",
    author_email="xikanfeng2@gmail.com",
    description="I-Impute: a coherent strategy to impute singlecell RNA sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xikanfeng2/I-Impute",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
