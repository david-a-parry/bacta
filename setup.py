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
          'parse_vcf>=0.2.1', #for bin/var_vs_contam
      ],
    scripts = ["bin/bacta", "bin/contam_reads_to_fastq", "bin/summary_to_bed",
               "bin/var_vs_contam",
               "bin/test_bacta_clipping_thresholds",
               "bin/test_bacta_expect_thresholds",
               "bin/remove_clipped_reads",
        ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        ],
)
