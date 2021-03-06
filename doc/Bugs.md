
# Kown bugs in Vulture


## FIXED IN 0.6.2

ID: 000005
Date: 20/02/2013
Severity: Normal
Priority: P3
OS: All
Assignee: ian.flinoft@york.ac.uk
Status: Needs further investigation
Resolution: Fixed. Conceptual error in scaling calculation.
Workaround: None
Summary: The direction vectors of the wave-vector and electric field
         in the planewave source may not always be rendered correctly in
         the gnuplot-planewave.dat file.

ID: 000004
Date: 20/02/2013
Severity: Normal
Priority: P3
OS: All
Assignee: ian.flinoft@york.ac.uk
Status: Confirmed
Resolution: Logic error in surface.c. Moved
            isInternalSurfaceType[BT_UNDEFINED] = true
            into addInternalSurface.
Workaround: Add manually to mesh.gnp file
Summary: The file containing internal surfaces, gnuplot-surface.dat, is
         not added to the list of files to splot in mesh.gnp.

## FIXED IN 0.6.1


ID: 000003
Date: 22/11/2012
Severity: Normal
Priority: P3
OS: All
Assignee: ian.flinoft@york.ac.uk
Status: Confirmed
Resolution: Fixed by extending ibbox by one cell in polarisation direction
            when rendering arrow.
Workaround: None
Summary: Incorrect rendering of current sources by gvulture.
         If size of bbox in source polarisation direction is zero
         arrow algorithm fails. Make it one cell long - half below, half above?

## FIXED IN 0.6.0

ID: 000001
Date: 19/11/2012
Severity: Normal
Priority: P2
OS: All
Assignee: ian.flinoft@york.ac.uk
Status: Fixed 20/11/2012
Resolution: Wrong times sent to incident field function and incorrect
            array indices in H field updates. 
Workaround: None
Summary: Getting strong reflection from TF/SF boundary.

ID: 000002
Date: 21/11/2012
Severity: Normal
Priority: P2
OS: All
Assignee: ian.flinoft@york.ac.uk
Status: Fixed 21/11/2012
Resolution: Wrong scaling function used to apply magnetic current density
            sources to grid.
Workaround: Scale amplitude in mesh file.
Summary: Magnetic current sources appear to have incorrect amplitude
         by factor equal to mesh size.
