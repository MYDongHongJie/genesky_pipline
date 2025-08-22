import sys
import json

with open(sys.argv[1]) as f:
    data = json.load(f)

summary = data["summary"]

before = summary["before_filtering"]
after = summary["after_filtering"]

metrics = [
    ("raw_reads", before["total_reads"]),
    ("raw_bases", before["total_bases"]),
    ("raw_q20", before["q20_rate"]),
    ("raw_q30", before["q30_rate"]),
    ("raw_gc", before["gc_content"]),
    ("clean_reads", after["total_reads"]),
    ("clean_bases", after["total_bases"]),
    ("clean_q20", after["q20_rate"]),
    ("clean_q30", after["q30_rate"]),
    ("clean_gc", after["gc_content"])
]

with open(sys.argv[2], "w") as out:
    for name, value in metrics:
        out.write(f"{name}\t{value}\n")
