# Protein-Ligand Movie Script


python
PDB_ID = "1hsg"
cmd.fetch(PDB_ID)
python end

hide everything
bg_color white

# Simple selections
select ligand, hetatm and not solvent
select protein_chains, polymer

# Visualization
show surface, protein_chains
set transparency, 0.5, protein_chains
color gray80, protein_chains

show sticks, ligand
color magenta, ligand

select binding_site, byres (ligand around 5) and polymer
show surface, binding_site
color yellow, binding_site
set transparency, 0.3, binding_site

# Movie setup
orient all
zoom all, 1.2
mstop
mclear
mset 1x360  

# Keyframes
frame 1
mview store, 1

frame 60
turn y, 180
mview store, 60

frame 120
turn y, 180
mview store, 120

frame 180
zoom ligand, 6
turn x, -20
mview store, 180

# Stay zoomed in from frame 180 to 360
frame 360
mview store, 360

# Play
set movie_loop, 1
mplay
