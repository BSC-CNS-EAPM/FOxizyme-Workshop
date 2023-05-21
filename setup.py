import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="foxizyme",
    version="0.0.1",
    author="EAPM",
    author_email="eapm@bsc.es",
    description="A module for setting up bioprospecting protocols for the FOxizymes workshop",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
    ],
)
