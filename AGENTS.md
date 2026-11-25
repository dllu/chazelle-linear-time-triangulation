# Chazelle's Linear Time Triangulation of a Simple Polygon

This project aims to fully implement Bernard Chazelle's linear time triangulation of a simple polygon in Rust.

Deliverables:

- fully working, clean, modular, maintainable Rust code
- maintain a DEVIATIONS.md that documents known bugs in the original papers, counterexamples, and their workarounds
- maintain a diary of your progress thus far in the `diary` directory in the form `yyyy-mm-dd-hh-mm-ss-title.md` in the local time zone
- you may consult an existing C++ implementation at https://github.com/Jayyhk/chazelle-triangulation. However, you must base your own implementation directly on the paper rather than merely transcribing the C++ code. That github has some transcribed copies of the original papers. You can also do your own search and cross reference to see if the transcriptions are accurate.
- maintain a minimal README.md that documents the current state and how to run it etc. Be mindful to avoid accumulating tons and tons of stale results or excessive details in the README.md by refactoring it and cleaning it up frequently.
- experimental evidence showing the O(n) performance compared to O(n log n) methods, even if there is a large constant factor. Add it as a chart in the readme.
