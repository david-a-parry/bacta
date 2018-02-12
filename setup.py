from distutils.core import setup
import glob

setup(
    name = "bacta",
    packages = ["bacta"],
    version = "0.0.1",
    description = "TO DO",
    author = "David A. Parry",
    author_email = "gantzgraf@github.com",
    url = "https://github.com/gantzgraf/bacta",
    license='GPLv3',
    install_requires=[
          'pysam>=0.12.0',
      ],
    scripts = ["bin/bacta", "bin/contam_reads_to_fastq", "bin/summary_to_bed",
               "bin/test_bacta_clipping_thresholds", 
                "bin/test_bacta_expect_thresholds",
        ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        ],
)
