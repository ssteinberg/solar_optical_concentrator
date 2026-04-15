# Differentiable Rendering
## Run Instructions
1. Install Python3
2. Install the requirements listed in the `requirements.txt`
   - If you want to use double precision, or any of the non-standard variants, you will need to compile your own version of Mitsuba3 locally rather than use an installed option (instructions may be found [here](https://mitsuba.readthedocs.io/en/stable/src/developer_guide/compiling.html))
3. Set the configuration variables at the top of the file as desired (pay special mind to the variant as you will have to select one your system supports)
4. Run the Jupyter Notebook as usual

## Notebooks
`single_mirror_spline` and `double_mirror_spline` are the main files to run, as they have the most robust set-ups (`single_mirror` and `double_mirror` are mostly kept for historical purposes).  
`single_mirror_spline` accounts for the case where a single, parabolic mirror is focussing onto a point (a line due to the 3D extrusion).  
`double_mirror_spline` accounts for the case of a two mirror concentrator, attempting to replicate the results of the paper.  
Please note that the `double_mirror_spline` was worked on last so it has some improvements not found in earlier notebooks.
