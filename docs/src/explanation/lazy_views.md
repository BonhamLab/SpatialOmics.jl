# Lazy views and spatial filtering

## The problem: filter-on-access, not filter-on-define

A 50-million-transcript Xenium dataset takes several seconds to load from disk.
Filtering it to a 500 µm × 500 µm region of interest requires reading only the
points within that bounding box — but only once you know which region you want.
If you filter eagerly at definition time, you pay the full read cost immediately,
lose the ability to redefine the region without re-reading, and produce a copy
that diverges from the on-disk source.

SpatialOmics follows the `SubArray` design from Julia's base library: define
the view cheaply, pay the cost only when the data is actually needed.

## SpatialExtent and SpatialROI

[`SpatialExtent`](@ref) is a typed axis-aligned bounding box. It carries `coord_system`
in addition to `xmin/xmax/ymin/ymax`, so the package can catch the common
mistake of filtering in the wrong coordinate space — a mismatch raises an error
rather than producing silently wrong results.

[`SpatialROI`](@ref) wraps any GeoInterface-compatible polygon with a cached bounding
box. The cache is used as a pre-filter: only points or shapes whose bounding box
intersects the ROI extent are tested for exact polygon containment, keeping the
cost proportional to the density of candidates rather than the total dataset
size.

Both types are accepted by `view`. They can also be converted to [`SpatialShapes`](@ref)
for visualisation — `SpatialShapes(ext)` produces a rectangular polygon.

## SpatialElementView and SpatialDatasetView

`view(el, roi)` returns a [`SpatialElementView`](@ref) — a struct holding a
reference to the parent element and the ROI. No data is read, no arrays are
allocated. The element's accessors — [`coords`](@ref), [`geometries`](@ref), [`feature_ids`](@ref),
[`instance_id`](@ref), [`count_per_instance`](@ref) — are all defined on [`SpatialElementView`](@ref)
and apply the filter on each call.

`view(ds, roi)` returns a [`SpatialDatasetView`](@ref) that applies the same ROI
to every element in the dataset. Accessing a specific element via
[`points`](@ref) or [`images`](@ref) returns a [`SpatialElementView`](@ref)
for that element.

This is the preferred way to build multi-layer plots: define the region once,
then pass the view to each plot verb independently. Each verb applies the filter
at materialisation time:

```julia
roi = view(ds, SpatialExtent(1000.0, 2000.0, 500.0, 1500.0; coord_system="global"))

fig = Figure(size=(600, 600))
ax  = Axis(fig[1, 1]; aspect=DataAspect(), yreversed=true)
heatmap!(ax, images(roi, "morphology"))        # crops image to extent
poly!(ax,    shapes(roi, "cell_boundaries"))   # filters shapes by extent
scatter!(ax, points(roi, "transcripts"))       # filters points by extent
fig
```

## When to use collect

`collect(view(el, roi))` materialises the view into a concrete element — a new
[`SpatialPoints`](@ref) or [`SpatialShapes`](@ref) with only the filtered rows. Use `collect`
when:

- You need a standalone element for downstream computation (e.g., passing to
  [`analyze`](@ref)).
- You are going to access the filtered data many times and want to avoid
  recomputing the mask on each call.

`collect` is an explicit opt-in. Keeping the view avoids the copy cost when
you only need to plot or inspect the region once.

## Overlap semantics for shapes

Point containment is unambiguous — a point is either inside a region or not.
Shape containment admits two interpretations:

- `:any` (default) — include shapes whose bounding box intersects the ROI.
  This is a fast approximation: some included shapes may extend outside the ROI.
- `:full` — include only shapes whose bounding box lies entirely within the ROI.
  More conservative; use when you need all included shapes to be completely visible.

The `overlap` keyword is passed to `view(el, roi; overlap=:any)`. It has no
effect on point filtering.
