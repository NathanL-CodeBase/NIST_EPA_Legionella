# Report review log, 2026-09-10

Source: `Legionella Relevant Aerosol Emissions from Showers to Adjacent Rooms.docx`
Output: `Legionella Relevant Aerosol Emissions from Showers to Adjacent Rooms_2026-09-10.docx`

The original is untouched. All edits in the new file are Word tracked changes authored as
"Claude (review)" so they can be accepted or rejected individually. All open questions are
Word comments from the same author. The document's 28 images, 10 embedded fonts, and all
equation objects were carried across without re-encoding; every part except
`document.xml`, `comments.xml`, `commentsExtended.xml`, `commentsIds.xml`,
`commentsExtensible.xml` and `people.xml` is byte-identical to the original.

Seven tracked edits and 15 comments. The two existing comments from Dustin Poppendieck are
preserved.

---

## Tracked changes

| # | Section | Change | Basis |
|---|---|---|---|
| 1 | Water Flow Rates | "The density of the was assumed to be 1 L/Kg" to "The density of the water was assumed to be 1 kg/L" | Missing noun, and the unit was inverted: 1 kg/L is the density of water, 1 L/kg is its specific volume |
| 2 | Experimental Procedure | Shower protocol time 02:00 to 03:00 | Confirmed by you 2026-09-10. `moduair_event_peak_times.py` gates ON transitions to ±5 min of 03:00 and 15:00, and the predefined event exclusions in `event_manager.py` are timestamped 15:00 |
| 3 | Experimental Procedure | Shower protocol time 14:00 to 15:00 | Same |
| 4 | Aerosol Monitoring | "23 bins sizes" to "24 bin sizes" | The Alphasense OPC-N3 reports 24 bins across 0.35 µm to 40 µm |
| 5 | Aerosol Monitoring | "only the first eight bins' sizes were used" to "only the first 12 bin sizes were used" | Confirmed by you 2026-09-10. The pipeline analyzes bins 0 to 11. The original also listed nine ranges under the word "eight" |
| 6 | Aerosol Monitoring | Bin list extended from 5.2 µm to 10.0 µm, adding 5.2 to 6.5, 6.5 to 8.0, and 8.0 to 10.0 µm | Follows from edit 5. Inserted runs reuse the document's existing Symbol-font pattern for µ so the styling matches |
| 7 | Table 6, emission parameters | Bedroom control volume 35.9 m³ to 36.1 m³ | Confirmed by you 2026-09-10. `BEDROOM_VOLUME_M3 = 36.1` in `src/particle_calculations.py`, annotated 36.10859771 m³ from CAD. 54.6 m³ was left alone and carries a comment instead |

---

## Comments raised

| # | Anchored at | Question |
|---|---|---|
| 1 | Fig. 1 caption | Caption numbering is mixed: 19 captions are SEQ fields, 17 are hard-typed. Full audit below |
| 2 | Aerosol Loss Rate, Section cross-reference | The REF field in "Tests performed" points at bookmark `_Ref237492571`, which sits on the Table 5 footnote rather than a heading |
| 3 | Air Change, λ values | 0.64 h⁻¹ ± 48 h⁻¹ and 0.93 h⁻¹ ± 89 h⁻¹ look like lost decimal points |
| 4 | Aerosol Loss Rate | βloss,entry is given as 0.93 h⁻¹ ± 0.89 h⁻¹, identical to λentry two paragraphs earlier |
| 5 | Fig. 9 caption | Onset slope 1.25 in the caption against 1.29 in the body text |
| 6 | Uniform Concentration Corrections | 73 experiments in one paragraph, 77 in the next and in the Fig. 10 caption |
| 7 | Aerosol Monitoring | Table 6 penetration factors still stops at 4.0 µm although the text now says 12 bins to 10.0 µm |
| 8 | Table 6, emission parameters | 54.6 m³ has no recorded derivation |
| 9 | Table 6, emission parameters | The eight emission variants, C_adjusted room, and both control volumes are implemented in no script in the repository |
| 10 | Equation 11 | Multiplying a mean of n rates by a single Δt yields a particle count, not a rate |
| 11 | Air Change exclusion criteria | Report says R² < 0.75, code applies 0.65, code docstrings say 0.75, and comment 82 in this document asks to raise it to 0.75 |
| 12 | Penetration Factor | Report defines p over a single 6 h to 1 h pre-shower window; the superseded code used paired 6 h windows before and after |
| 13 | Tracer Gas Injection System | Injection ran 4 min then 6 min from 2026-01-22; this section describes only the later protocol |
| 14 | Aerosol Monitoring | Every µm is a Symbol-font run, which extracts as "mm" in PDF text, copy and paste, and screen readers |
| 15 | Tests performed | Flooring changed in both rooms on 2026-05-13 and 2026-05-21, but only those two events are excluded, and βloss is a surface deposition term |

---

## Caption numbering audit

This is the one structural problem in the document, and it cannot be fixed by editing
numbers. 19 captions use Word `SEQ` fields and renumber on F9. 17 carry hard-typed numbers
and never will. The two sets have drifted apart, which is why three captions currently read
"Fig. 10" and three read "Table 6".

**Hard-typed, will not renumber:** Fig. 1, Fig. 2, Fig. 3, Fig. 4 (Experimental procedure),
Fig. 5, Fig. 6, Fig. 7, Fig. 8, Fig. 9, Fig. 10 (Cbed1 to Croom ratio), Fig. 15 (appendix,
bin 5 mixing data), Table 1, Table 2, Table 3, Table 5, Table 6 (elevation zones),
Equation 1.

**SEQ fields, renumber automatically:** Fig. 4 (Experimental temperature), Fig. 4
(Experimental relative humidity), Fig. 10 (air change rate), Fig. 10 (aerosol loss rate),
Fig. 11 through Fig. 14, Table 6 (penetration factors), Table 6 (emission parameters),
Equation 2, Equation 3 (twice), Equation 5, Equation 7 through Equation 11.

**Fix:** rebuild each hard-typed caption with References > Insert Caption, then select all
and press F9, then refresh the List of Tables and List of Figures. Do this before acting on
any figure or table cross-reference, because the numbers will move.

Also outstanding in the same area: Equation 4 and Equation 6 have no captioned paragraph.
Equation 6 exists in the Air Change section but its paragraph does not carry the Caption
style, so it is invisible to the sequence. Equation 4 appears to have been deleted.

---

## Deliberately not changed

**Every "mm" in the document.** These are not unit errors. Each µ is a Symbol-font run
containing the letter `m`, which renders correctly as µ in Word. There are 54 such runs.
Converting them to Unicode U+00B5 is worth doing for PDF text extraction and accessibility,
but it touches text that is currently correct on screen, so it is raised as comment 14
rather than done unasked.

**Figure, table, equation, and section numbers.** See the audit above. The "Section 3.4"
and "Section 3.4.1" cross-references are REF fields whose displayed values are stale caches,
not typos; they resolve on F9. Only the bookmark target in comment 2 is genuinely wrong.

**Every number listed in the comments.** Per your instruction on 2026-09-10, unresolved
values are raised as comments rather than corrected, so you can check them against the
analysis outputs first.

**The empty sections.** Abstract, Executive Summary, Acknowledgments, and Results are left
as they stand.

**The two Fig. 4 placeholder captions** for experimental temperature and relative humidity.
`scripts/hobo_onset_decay_figures.py` now produces exactly these as
`output/plots/hobo/hobo_temp_pre_post.html` and `hobo_rh_pre_post.html`, showing one pooled
pre point and one pooled post point per event across all five HOBO sensors with pooled
standard deviation whiskers, plus a paired t-test summary. You said you will insert them
when joining everything together.
