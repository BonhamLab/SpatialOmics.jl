_reserved_experiment_keys = (
    "slides"
)

_reserved_slide_keys = (
    "fovs"
)

reserved_fov_keys = (
    "name",
    "x",
    "y"
)


function load_experiment(name, exppath)

    ex = Experiment(name, exppath)
    slides_path = joinpath(exppath, name, "slides")

    for (sname, slide) in slides
        slide_path = joinpath(exppath, "slides", sname)
        for (i, (fname, fov)) in enumerate(slide["fovs"])
            f = FOV(i)
            f.properties["name"] = fname
            f.properties["x"] = get(fov, "x", NaN)
            f.properties["y"] = get(fov, "y", NaN)
            for (key, value) in pairs(fov)
                if key ∉ reserved_fov_keys
                    f.properties[key] = value
                end
            end
            push!()
        end
        

    end
    return ex
end


