#!/bin/bash

# ==============================================================================
# calc_electrostatics.sh
# 
# Description:
# Calculates the time-averaged 3D electrostatic potential of the α7-nAChR 
# using the Particle Mesh Ewald (PME) method via VMD's pmepot plugin.
# 
# The script centers the trajectory on the channel pore (VAL 251), computes 
# the potential on a 3D grid, and saves the output as a Data Explorer (.dx) file.
# It then uses awk to translate and time-average the .dx volumetric data into a 
# clean 4-column text format (X, Y, Z, Potential).
# ==============================================================================

# 1. Default Parameters
GRID=200
EWALD_FACTOR=0.25

# 2. Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --top) TOP="$2"; shift ;;
        --traj) TRAJ="$2"; shift ;;
        --out_prefix) OUT_PREFIX="$2"; shift ;;
        --grid) GRID="$2"; shift ;;
        --ewald) EWALD_FACTOR="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Validate required inputs
if [ -z "$TOP" ] || [ -z "$TRAJ" ] || [ -z "$OUT_PREFIX" ]; then
    echo "Usage: ./calc_electrostatics.sh --top <topology> --traj <trajectory> --out_prefix <prefix>"
    exit 1
fi

DX_FILE="${OUT_PREFIX}_pmepot.dx"
RAW_FILE="${OUT_PREFIX}_pmepot.rawdat"

echo "Calculating PME Electrostatics for ${TRAJ}..."

# 3. Generate the VMD Tcl script dynamically
cat << EOF > pmepot_temp.vmdin
package require pbctools
package require pmepot

# Load topology and trajectory
mol new $TOP
mol addfile $TRAJ waitfor all 0

# Center and wrap the system around the pore (VAL 251)
pbc wrap -center com -centersel "resname VAL and resid 251" -compound residue -all

set numframes [molinfo top get numframes]
puts "Aligning \$numframes frames..."

for {set i 0} {\$i < \$numframes} {incr i} {
    animate goto \$i
    set all [atomselect top "all"]
    set sel [atomselect top "resname VAL and resid 251"]
    set com [measure center \$sel weight mass] 
    
    # Align center of mass to origin
    \$all moveby [vecinvert \$com] 
    \$all delete
    \$sel delete
}

# Run the PMEpot plugin on all frames
puts "Running pmepot plugin..."
pmepot -frames all -grid ${GRID} -dxfile ${DX_FILE} -ewaldfactor ${EWALD_FACTOR}

exit
EOF

# 4. Execute VMD in text mode
vmd -dispdev text -e pmepot_temp.vmdin
rm pmepot_temp.vmdin

# 5. Process the output .dx file
# The .dx file contains the volumetric grid data. This awk script extracts the grid 
# boundaries, parses the data block, and computes the spatial coordinates (X, Y, Z) 
# and the time-averaged potential for each grid point.
echo "Processing .dx file and extracting time-averaged coordinates..."

awk '
BEGIN { ind = 0 }
# Extract grid dimensions
(NR==2) { Nx=$6; Ny=$7; Nz=$8 }
# Extract origin coordinates
(NR==3) { ox=$2; oy=$3; oz=$4 }
# Extract delta steps along X, Y, and Z
(NR==4) { dx=$2 }
(NR==5) { dy=$3 }
(NR==6) { dz=$4 }
# Parse the data block (Starts after line 8)
(NR>8) {
    for(c=1; c<=NF; c++) {
        ind++; 
        i = int(ind/Ny/Nz) % Nx; 
        j = int(ind/Nz) % Ny;
        k = (ind/1) % Nz; 
        sum[i,j,k] += $c;
        count[i,j,k]++
    }
}
# Calculate averages and output standard X Y Z Potential format
END {
    for (i=0; i<Nx; i++) {
        for (j=0; j<Ny; j++) {
            for (k=0; k<Nz; k++) {
                if (count[i,j,k] > 0) {
                    print ox + dx*i, oy + dy*j, oz + dz*k, sum[i,j,k]/count[i,j,k]
                }
            }
        }
    }
}' < ${DX_FILE} > ${RAW_FILE}

# 6. Cleanup large temporary files
echo "Cleaning up temporary files..."
rm ${DX_FILE}

echo "Electrostatic calculation complete. Data saved to ${RAW_FILE}"