scale_colour_viridis <- function (..., type = "seq", palette = 1) 
{
    discrete_scale("colour", "viridis", viridis, 
                   ...)
}

scale_fill_viridis <- function (..., type = "seq", palette = 1) 
{
    discrete_scale("fill", "viridis", viridis, 
                   ...)
}