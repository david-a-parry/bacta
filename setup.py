from distutils.core import setup
import glob

setup(
    name = "bacta",
    packages = ["bacta"],
    version = "0.0.1",
    description = "TO DO",
    author = "David A. Parry",
    author_email = "gantzgraf@github.com",
    url = "",
    scripts = glob.glob( "bin/*.py"),
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        ],
)