# Assign transcripts and summarise expression

```@setup expression-summaries
using CairoMakie
using Markdown
CairoMakie.activate!(type="svg")
set_theme!(Theme(
    fontsize=15,
    Figure=(; backgroundcolor=:white),
    Axis=(; xgridvisible=false, ygridvisible=false),
))
```

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

The spatial join and its matrix summary are two views of the same result. The
unmatched transcript remains visible spatially but contributes to no matrix
row.

```@eval expression-summaries
figure = Figure(size=(850, 370))
spatial_axis = Axis(figure[1, 1]; aspect=DataAspect(), title="Spatial membership",
                    xlabel="x (µm)", ylabel="y (µm)")
poly!(spatial_axis, geometries(cells); color=(:lightsteelblue, 0.4),
      strokecolor=:steelblue, strokewidth=2)
gene_palette = [:darkorange, :seagreen]
for gene_id in eachindex(features(transcripts))
    mask = feature_ids(transcripts) .== gene_id
    scatter!(spatial_axis, coords(transcripts)[mask]; color=gene_palette[gene_id],
             markersize=15, label=features(transcripts)[gene_id])
end
text!(spatial_axis, 2, 3.45; text="cell 101", align=(:center, :center))
text!(spatial_axis, 8, 3.45; text="cell 102", align=(:center, :center))
text!(spatial_axis, 21, 20.5; text="unmatched", align=(:right, :bottom), color=:gray35)
axislegend(spatial_axis; position=:lt, framevisible=false)
xlims!(spatial_axis, -1, 22); ylims!(spatial_axis, -1, 22)
matrix_axis = Axis(figure[1, 2]; title="Cell-by-gene counts",
                   xticks=(1:2, var_names(expression)),
                   yticks=(1:2, string.(source_ids(expression))),
                   xlabel="gene", ylabel="cell instance ID", yreversed=true)
matrix_values = Matrix(expression[source_ids(expression), var_names(expression)])
heatmap!(matrix_axis, matrix_values; colormap=:Blues, colorrange=(0, maximum(matrix_values)))
for row in axes(matrix_values, 1), column in axes(matrix_values, 2)
    text!(matrix_axis, column, row; text=string(matrix_values[row, column]),
          align=(:center, :center), color=:black)
end
save("expression-summaries.svg", figure)
Markdown.parse("![Spatial memberships and cell-by-gene counts](expression-summaries.svg)")
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
