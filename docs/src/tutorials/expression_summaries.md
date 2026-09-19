# Assign transcripts and summarise expression

Spatial assignment is a many-to-many spatial join. Unmatched transcripts are
absent from a membership relation, while a point contained by overlapping cell
objects can appear more than once. This makes the cardinality explicit instead
of silently choosing an object or inventing a sentinel ID.

## Construct transcripts and cells

```@example expression-summaries
using SpatialOmics

square(xmin, xmax, ymin, ymax) = Polygon(Point2f[
    (xmin, ymin), (xmax, ymin), (xmax, ymax),
    (xmin, ymax), (xmin, ymin),
])

transcripts = SpatialPoints(
    [Point2f(1, 1), Point2f(2, 2), Point2f(7, 1), Point2f(20, 20)];
    feature_id=Int32[1, 2, 1, 2],
    feature_codebook=["Actb", "Gapdh"],
    coord_system="global_um",
)
cells = SpatialShapes(
    [square(0, 4, 0, 4), square(6, 10, 0, 4)];
    instance_id=Int32[101, 102],
    coord_system="global_um",
)
nothing
```

## Inspect transcript-to-cell membership

For point membership, `source_ids` are one-based rows of the point collection;
`destination_ids` are cell `instance_id` values. The fourth transcript is
outside both cells and therefore has no relation row.

```@example expression-summaries
membership = analyze(Membership(), transcripts, cells)
(
    transcript_rows=source_ids(membership),
    cell_ids=destination_ids(membership),
    matched_fraction=length(unique(source_ids(membership))) / length(transcripts),
)
```

Do not compute coverage as `nobs(membership) / length(transcripts)` when cell
objects may overlap: one transcript can contribute multiple membership rows.
Multipolygon components belonging to the same cell are deduplicated.

## Build a cell-by-gene matrix

`Expression` uses the same indexed point-in-polygon join and accumulates a
cell-by-gene count matrix. Rows are selected by cell instance ID and columns by
feature name.

```@example expression-summaries
expression = analyze(Expression(), transcripts, cells)
(
    cell_ids=source_ids(expression),
    genes=var_names(expression),
    cell_101=expression[101, :],
    actb=expression[:, "Actb"],
)
```

The result is intentionally a lightweight relation rather than a full
single-cell analysis object. Use `annotate` for row metadata, then hand the
matrix and metadata explicitly to clustering, normalization, or dimensionality
reduction packages.

```@example expression-summaries
count_matrix = expression[source_ids(expression), var_names(expression)]
size(count_matrix)
```

```@example expression-summaries
annotated = annotate(expression, ["left", "right"]; key=:region)
(nobs(annotated), nvar(annotated), obs_names(annotated))
```

`obs_names` falls back to the cell instance IDs because this annotation used
the key `:region`; use `key=:name` when string row indexing is desired.
