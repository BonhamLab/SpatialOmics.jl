module NearestNeighborsExt

using NearestNeighbors
using Tables
using SpatialOmics

function SpatialOmics.analyze(tok::KNN, rel::SpatialRelation{Expression})
    nobs(rel) == 0 && error("Cannot compute KNN on empty relation")

    obs = Tables.columntable(rel.obs)
    if hasproperty(obs, :x) && hasproperty(obs, :y)
        xs, ys = Float32.(obs.x), Float32.(obs.y)
    elseif hasproperty(obs, :centroid_x) && hasproperty(obs, :centroid_y)
        xs, ys = Float32.(obs.centroid_x), Float32.(obs.centroid_y)
    else
        error("KNN requires obs columns :x/:y or :centroid_x/:centroid_y. " *
              "Add centroid coordinates via annotate(rel, coords; key=:x) first.")
    end

    n     = length(xs)
    k_eff = min(tok.k, n - 1)
    k_eff < 1 && error("KNN requires at least 2 observations (got $n)")

    pts = Matrix{Float32}(undef, 2, n)
    pts[1, :] = xs
    pts[2, :] = ys

    tree            = KDTree(pts)
    idx, raw_dists  = knn(tree, pts, k_eff + 1, true)   # +1 includes self

    dst_ids = zeros(Int32, n * k_eff)
    weights = zeros(Float32, n, k_eff)

    for i in 1:n
        nbrs  = idx[i]
        ndist = raw_dists[i]
        skip  = raw_dists[i][1] ≈ 0f0 ? 1 : 0   # skip self if at distance 0
        for j in 1:k_eff
            src_j = j + skip
            src_j > length(nbrs) && break
            weights[i, j]            = Float32(ndist[src_j])
            dst_ids[(i-1)*k_eff + j] = rel.src_ids[nbrs[src_j]]
        end
    end

    SpatialRelation(KNN(k_eff), rel.src, rel.src,
                    rel.src_ids, dst_ids, weights; obs=rel.obs)
end

end
