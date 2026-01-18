#!/usr/bin/env nu

$env.GITROOT = (git rev-parse --show-toplevel | str trim)
let neb_base = ($env.GITROOT | path join "eonRuns" "results" "03_neb")

# Get all system directories and process them in parallel
ls $neb_base | where type == dir | par-each { |sys|
    let sys_name = ($sys.name | path basename)
    
    ["mmf", "cineb"] | each { |method|
        let work_dir = $"($sys.name)/($method)"
        if not ($work_dir | path exists) { return }

        # Match the "other" method for the overlay
        let cfg = if $method == "mmf" {
            { other: "cineb", label: "CINEB", out: "cineb_on_roneb" }
        } else {
            { other: "mmf",   label: "RONEB", out: "roneb_on_cineb" }
        }

        print $"Plotting ($sys_name) with ($method)..."
        
        cd $work_dir
        python -m rgpycrumbs.cli eon plt_neb --con-file neb.con --plot-structures "crit_points" --facecolor "floralwhite" --plot-type landscape --figsize 5.37 5.37 --dpi 300 --zoom-ratio 0.25 --fontsize-base 12 --ase-rotation 0x,0y,0z --additional-con $"../($cfg.other)/sp.con" $cfg.label --highlight-last --title "" --show-legend --ira-kmax 14 -o $"../../../04_plots/($sys_name)/($method)/($cfg.out).png"
        cd -
    }
}

print "All reactions plotted successfully."
