# Vulture To Do List

## Priority Items

### Generalise tests to use Octave or MATLAB. 

This will require turning all Octave sctipts into functions- same as AEG Mesher.

### Find way to reduce amount of static validation data in tests

### Implement time domain utilities as MATLAB functions

Implement xtime, xplane, mxplane, xtransall, xfreq, xfplane and mxfplane as 
Octave/MATALB functions which can also output in physical units and possibly 
interpolate E/H onto nodes and integer time steps. Issues are edges of sampled 
space/time.

### Add mesh origin directive

    MO <real: xo> <real: yo> <real: zo>
    
Mesh origin directive for use with MS in case of cubic/cuboid mesh so physical 
coordinates are still correct relative to input mesh. Default to (0,0,0).

### Add basic thin wire model and delta gap voltage source

See section wires below.

### Optimise mesh parsing

Mesh parsing may be using too much RAM when each cell/face is in separate MB/TB. 
Can this be optimised?

### Completely rethink observers and interaction with mesher

See below.


## Build and packaging infrastructure

### Support older cmakes

### Improve support for different compilers


## Documentation

### Write script to generate figures in manual from test-suite results

### Add Implementation Manual once stablised

### Update Implementation Manual on plane-waves, interations of media/boundaries


## Testsuite

### Parser

 * `parser_eof` for EOF issue.

 * Undefined tags.
 
 * Missing required directives.
 
 * Duplicated singleton directives.

### Sources

 * softsrc_ix/iy/iz, soft electric current source IX/IY/IZ
 * softsrc_ex/ey/ez, soft electric field source EX/EY/EZ
 * softsrc_imx/imy/imz, soft magnetic current source IX/IY/IZ
 * softsrc_hx/hy/hz, soft magnetic field source  HX/HY/HZ
 * hardsrc_hx/hy/hz, hard magnetic field source with integrated waveform.
 * lumpedload_vrx/vry/vrz, zero voltage - just load.
 * hard_src_jmsxy, aperture radiation by impressed surface magnetic current (Balanis)
 * hard_src_jmsyx
 * hard_src_jmsxz
 * hard_src_jmszx
 * hard_src_jmsyz
 * hard_src_jmszy

### External BC

 * pml_pecwall_xxx, PML cutting PEC - no transmission
 * parplate1D_periodic_xxx, empty waveguide using PERIODIC BC.
 * microstrip_reflect_xxx, PML terminating inhomogeneous media:

### Graded mesh

 * parplate1D_gradedslab_xxx, simple dielectric slab in waveguide using graded mesh 
   with two bisections. Three directions.
 * aperture_x/y/z, grading of square aperture edge c.f polarisability `m = alpha_m * H_sc`.

### Internal surfaces

 * openbox_pec, hertzian dipole in PEC box with aperture - 
 * openbox_se_pec, ILCM model

### Planewaves

 * planewave_empty_oblique, random direction
 * planewave_aperture, radiated field from aperture, cf magnetic dipole. mag + phase

### AVERAGED_MEDIA

Test average media support for dielectric block in TEM waveguide. Compare BC with and 
without averaged medium support.

### Improve data comparision in pydifffd

Cope with

 1. Zero crossings.
 2. Noise floor.
 3. Deviations above Nyquist or lambda/10 frequency.

Reactive in tests once reliable.

### Add timing folder with timing tests

Empty grid, hertzian dipole source to excitate al polarisations, single point ASCCI observer, Mur ABC.

    function( vulture_timing_test , TESTNAME , EXTENT ) 

      configure_file( ${VULTURE_SOURCE_DIR}/timing/timing.mesh.cmake ${VULTURE_BINARY_DIR}/timing/${TESTNAME}.mesh )

      add_test( NAME ${TESTNAME}_vulture COMMAND ${VULTURE_BINARY_DIR}/src/vulture -v ${TESTNAME}.mesh )

      if( PYTHONINTERP_FOUND )
        add_test( NAME ${TESTNAME}_extract COMMAND ${PYTHON_EXECUTABLE} ${VULTURE_SOURCE_DIR}/util/timing.py vulture.log results.dat )
      endif( PYTHONINTERP_FOUND )

    endfunction()

    vulture_timing_test( timing10 , 10 )

    ..........

Use function to create tests - configured source?

Python to extract data from log?

Get RAM size and limit max extent?

  timing/empty10
  timing/empty20
  timing/empty30
  timing/empty40
  timing/empty50
  timing/empty60
  timing/empty70
  timing/empty90
  timing/empty100
  timing/empty200
  timing/empty300
  timing/results.txt
 
     # Vulture: ?.?.?
     # Date: ??/??/????
     # Time: ??:??:??
     # Machine name: ??????
     # OS: ??????
     # Number of processors: ?
     # CPU: ?????
     # RAM: ?????
     # Compiler: ??????
     # Compile flags: ?????
     # Compile defines: ?????
     # Number of OpenMP threads: ?
     #
     # Model   Run   SPI    SPI/cell
     # size    Time  
     # [cells] [s]   [ms]   [ns]
     #------------------------------
     10        ?     ?      ?

     
## Interfacing

### Write falcon2vulture to translate a Falcon tdfd.in file into a Vulture mesh file

Python?

### Write hawk2vulture to translate a Hawk tlm3d.in file into a Vulture mesh file

Python?

### Write an AMELET wrapper


## Optimisation

### Pointers versus indices

Should sources/blocks/surfaces/wires hold pointers to their 
waveforms/medium/boundary/wires type rather than numbers which have to be 
searched for? Media need numbers for indexing but blocks don't really need them.

### Investigate run time performance cf UGR-FDTD on HIRF HWTC

Output flushing?

### Add dereferencing variables to TF/SF boundary updates

### Deallocate unrequired stuff at end of init


## External tools

### Extend processing code to allow multiple user defined excitations


## Parser

### Parser doesn't always find problems due not spotting extraneous junk at end of line.

### Default option treatment in parser not always good.

Should all be set in parseXX - not in init functions?


## Graphics

### gvulture: option to render observers as boxes rather than points

## Determine how to view data on a mesh using gmsh. 

Can gvulture be modified to translate impulse.dat output into op<n>.msh so it 
can be viewed using

    gmsh mesh.msh op<n>.msh?
    
### Still PW arrow scaling issue in gvulture with large aspect bbox

For example, 1D WGs. See if using min. box side length to determine arrow length is better. 

May still be duff is some cases - check aperture example and all pw test meshes.


## Grid

### Check mesh line are monotonically strictly increasing sets

### Refactoring: Make standard update equations into separate functions

    updateEfieldStandard( int flim[6][6] )
    updateHfieldStandard( int flim[6][6] )

    void updateGridEfield( void )
    {
      updateEfieldStandard( gfilim );
      return;
    }

    void updateGridHfield( void )
    {
      updateHfieldStandard( gfilim );
      return;
    }

    void setGrid( int flim[6][6] , real initialE , real initialH )
    {
      .....
    }

    void clearGrid( void )
    {
      setGrid( gfglim , INITIAL_FIELD_VALUE , INITIAL_FIELD_VALUE );
      return;
    }


## Media

### Fix average media support that was hosed by Debye related updates

### Is averaged media support correct with regard to PEC media

Should we average PEC or apply it directly after averaging. At the moment it takes 
part in the averaging using a conductivity of 1e8. PEC material type is in fact 
bypassed when average media is on.

### What about lumped resistive loads

### What about debye media?

### Tests need to deal with USE_AVERAGED_MEDIA

Will generate different answers depending on whether it is ON/OFF.

### Extend  Debye media into PML.

When init blocks if surface in on PML boundary extend it to gobox. But edge to 
extend into edge and corner regions of overlapping PML faces so not trivial.
 
### Voxel object loader

vox file

    <s: dataFile>
    <i: bytes>
    <i:nx> <i:ny> <i:nz>
    <r:dx> <r:dy> <r:dz>
    <r:xref> <r:yref> <r:zref>
    <i: key1> <r: eps_r1> <r: sigma1> <r: mu_r1>
    .....
    <i: keyN> <r: eps_rN> <r: sigmaN> <r: mu_rN>    

    <s: dataFile> is little little binary file with 
    <i: nx > * <i: nx > * <i: nx > unsigned integers of size  <i: byte>
 
stored in row major format

    for( i = 0 ; i < nx ; i++ ) 
      for( j = 0 ; j < ny ; j++ ) 
        for( k = 0 ; k < nz ; k++ ) 
          fread(   voxel[i][j][k] , 1 )

    MT normal VOXEL "norman.vox"
    MB 20 20 30 30 40 40 normal theta phi

Located reference point in vox file at node given. Rotations applied according 
to two angles. Object truncated at mesh boundaries

Transform voxel object bbox onto mesh. Located new grid aligned bbox, trruncated at 
mesh boundaries. Some of this technology already exists in MTHR project files.

#### Method 1 (Done at preprocesssing stage)

Transform cells centres within grid aligned bounding box back to object and use 
nearest neighbour interpolation to decided material of cell. Add MB for cell. 
Could be done for both indexed and averaged media but will likely generate a 
vast number of MBs.

#### Method 2 (Done at init stage, indexed media only)

Transform all edge centres within grid aligned bounding box back to object and 
use nearest neighbour interpolation to decided material of edge. Add material to 
material list and apply index directly to grid

#### Method 3 (Done at init stage, averaged media only)

Transform particle swarm around each edge centre within grid aligned bounding 
box back to object and use "particle number weighted" interpolation to get 
material of parameters of edge. Apply material coefficients directly to grid

How to view object on mesh with method 2 and 3?


## Sources

### Source meshing strategy
         
In the input unstructured mesh: Sources are normal group of elements - any dimensionality.

### New per group meshing option: physicalType = 'SOURCE'

Sources mapped normally according to dimensionality and setting of type.
In mesh exporter catch physicalType = 'SOURCE' and export source AABB in EX:
    
    WF wf1 COMPACT_PULSE 
    EX <ilo> <ihi> <jlo> <jhi> <klo> <khi> <groupName> EZ wf1
    
If no source found output plane wave 3 cells inside CV:

    PW <ilo> <ihi> <jlo> <jhi> <klo> <khi> <groupName> wf1 90 0 90
 
### Port PW masks to use FaceMask utilities

### Lumped scattering parameters

### Dipole moment sources: PX, PY, PZ, MX, MY, MZ

### Check if sources applied on PEC:

With PEC implemented through media coefficients important not to apply source on PEC locations.
Add checks in initSource that not applying source to field with PEC coefficients and give warning.

### SWE sources

    SW <ilo> <ihi> <jlo> <jhi> <klo> <khi> <name > SWE <r: x0> <r: y0> <r: z0> <i: order> <r: Qn_real> <r: Qn_imag> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]
    SW <ilo> <ihi> <jlo> <jhi> <klo> <khi> <name > SWE <r: x0> <r: y0> <r: z0> <s: fileName> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]

### Port sources

    EX <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT COAXIAL     <r: innerRadius>/<r: Zc> <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]
    EX <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT CIRCULAR    <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]
    EX <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT RECTANGULAR <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]

 * Bounding box must be surface.
 * Outer port size taken from bounding box.
 8 Modes: "TEM" "TE:1,0", "TM:1,1",....
 8 Excite only E or E/H to give one-way wave?  

### Can PW auxiliary grid approach be used on uniform non-cubic grid?

### Advanced TF/SF boundary

Inverse FFT method TF/SF

 * run empty grid with source only (vulture -i fred.mesh) - how to excite?
 * collect and store frequency reponse of points on Huygens box.
 * run model with scatterer (vulture fred.mesh).
 * loaded frequency response of huygens points.
 * use running iDFT to replay time response at points as required by TF/SF boundary.

Tan's stuff?


## Observers

### Observer meshing strategy

In Vulture EX, PW, OP can only be defined for cuboid (possibly degenerate) 
shapes. In more complex cases may want to e.g define observation of current on 
curve surface. Related to issue of how to get outputs out and map back onto 
input mesh.

In the input unstructured mesh: Observers are normal group of elements - any dimensionality.

New per group meshing option: physicalType = 'OBSERVER'

Observers mapped as all nodes in group using float indices.

In mesh exporter catch physicalType = 'OBSERVER' and export:
    
    OT <groupName> <format> <domain> <quantity1> ... <quantityN>
    OP <ilo> <ihi> <jlo> <jhi> <klo> <khi> <groupName> <param1> ... <paramN>     

    If no observers found output far-field observer 3 cells inside CV

    FF .......

    One OT per group, could be many OP per group. 

Order OPs according to order in mesh so the output data can be mapped later. 
    
For meshing complex surfaces each node is output separately, e.g.

    OT wing HDF5 FREQ E H
    OP 21.6 23.1 45.3 21.9 23.1 45.4 wing   
    OP 21.6 23.1 46.3 21.8 23.2 46.6 wing
    .........................
    OP 21.6 23.1 47.3 21.3 23.4 47.6 wing  

But for compatibility we still allow AABB with strides:

    OT wing HDF5 FREQ E H
    OP 2 2 2 10 10 10 wing 5 5 5

where 5 5 5 are the number output points in the range 2 to 10. The default is 1 so
only 2 would be output. stride_i = ( xhi - xlo ) / ( num_i - 1 ). 

Notes:

1. Ideally outputs should be keyed to input mesh and physical units.
   Favours nodal observers as the node positions are directly available in
   the input mesh.

2. Nodal observers require two point averaging of E and four point averaging
   of H. This will incur some performance hit - however - this would have to be done
   in post-processing anyway so fields are at known locations.
     
3. For bounding box will get a different number of points than for cell based 
   binary observer. This has implications for AEG processing output format.

4. Temporal averaging of magnetic fields is harder - need to cache last values
   so time-average can be taken H^n = 0.5 * ( H^{n-1/2} + H^{n+1/2} ). Could consume
   a large amount of RAM for volume observers with small steps.

5. Is caching a useful strategy for improving performance. Each time domain 
   observer caches, say, 20 time steps of data and then flushes to disk. Temporal 
   averaging of H could be done on flush (flush to last but one time step). Need to 
   retain last value for next batch - circular buffers with head and tail? Does not 
   generalise to other observables like S which are combination of E/H. Need to 
   keep temporal averaging of H and caching of observable separate. Would make 
   timer prediction erratic if cache interval >~ timer update interval?

6. Portable binary output format: hdf5? AMELET? May as well be AMELET? BUT HDF5 dosen't 
   currently support concurrent write/read so cannot monitor time responses 
   from another program. See:

   * http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2010-June/003211.html
   * http://www.hdfgroup.org/hdf5-quest.html#grdwt     

   How to write data into h5 file as it is generated?

   subject IDs:

     0  electricField
     1  magneticField
     2  powerDensity
     3  planeWaveDecomposition
     4  current
     5  voltage
     6  power
     7  sParameter
     8  zParameter
     9  yParameter
     10 theveninVoltageGenerator
     11 nortonCurrentGenerator
     12 couplingCrossSection
     13 radarCrossSection

     data.h5
     |-- label/
     |   `-- predefinedOutputRequests
     |-- mesh/
     |   `-- $gmesh1
     |       `-- $sphere
     |           `-- group
     |               |-- $inside[@type=volume]
     |               `-- $skin[@type=face]
     |-- floatingType/
     |   `-- $e_field
     `-- outputRequest/
         `-- $outputRequest_group/
             `-- $or1[@subject=/label/predefinedOutputRequests
                      @subject_id=0
                      @object=/mesh/$gmesh1/$sphere/group/$inside
                      @output=/floatingType/$e_field]
 
7. Only implement observers that cannot easily be done in post-processing:

   * Z -> adds V/I observers and flags impedances as post-processing.
   * S -> addes E/H observers and flags S as post processing.

9. Should there be sepatate E and H observers? Would reduce resource consumption
   if only one required.

10. Observer types:
 
     OP <int: XLO> <int: XHI> <int: YLO> <int: YHI> <int: ZLO> <int: ZHI> ...
        <tag: name> <id: type> [<id: refWaveform>]
     OP <int: XLO> <int: XHI> <int: YLO> <int: YHI> <int: ZLO> <int: ZHI> ...
        <tag: name> <id: type> <int: XSTEP> <int: YSTEP> <int: ZSTEP>  [<id: refWaveform>] 
     FF <int: XLO> <int: XHI> <int: YLO> <int: YHI> <int: ZLO> <int: ZHI> ...
        <tag:name > <real: theta1> <real: theta2> <int: num_theta> <real: phi1> <real: phi2> <int: num_phi> [ <mask> [<id: refWaveform>] ]

Mesh element              | Domain    | Format |Physical quantity
:-------------------------|:---------:|:------:|:----------------
node,line,surface,volume  | TIME      | AEGBIN | EH
 node                     | TIME      | ASCII  | EH,E,H
node                      | FREQ      | ASCII  | EH,E,H,S,P
 node,line,surface,volume | TIME      | HDF5   | E,H
node,line,surface,volume  | FREQ      | HDF5   | E,H,S,P
surface                   | TIME,FREQ | HDF5   | I
line                      | TIME,FREQ | HDF5   | V
volume                    | FREQ      | HDF5   | P,RCS


     FIELD_TDOM_ASCII - Single node ASCII time domain output.
                        Fields output at node position at electric field times.
                        Spatial averaging of E/H and temporal averaging of H.
                        Field stored in ASCII file <name>.asc.
             
     FIELD_FDOM_ASCII - Single node ASCII frequency domain output.
                        Fields output at node position with correct phase relative to input waveform.
                        Spatial averaging of E/H and temporal averaging of H then running DFT.
                        Field stored in ASCII file <name>.asc.
                   
     FIELD_TDOM_IMPULSE - bbox is node/line/surface/volume with defined step.
                          Cell based binary time domain output - processing tool impulse.dat compatible.
                          Outputs fields at their respective locations in the cell at their repsective times.
                          Field stored in impulse.dat with metadata in excite.dat/process.dat.

     FIELD_TDOM_HDF - bbox is node/line/surface/volume.
                      Node based binary time domain output. 
                      Outputs field at bbox nodes at electric field times.
                      Spatial averaging of E/H and temporal averaging of H.
                      Data stored in HDF file output.h5 at location /floatingType/<name>

     POWER_FDOM_HDF - bbox is node/line/surface/volume with step.
                      Power density sigma <|E|^2> is determined at nodes.
                      Is sigma obtainable from alpha/beta/gamma?
                      Frequency domain output using running DFT.
                      Data stored in HDF file output.h5 at location /floatingType/<name>

     VOLTAGE_FDOM_HDF - bbox is line.
                        Voltage difference integrated along line.
                        Frequency domain output using running DFT.
                        Data stored in HDF file output.h5 at location /floatingType/<name>

     CURRENT_FDOM_HDF - bbox is surface.
                        Current flow through surface.
                        Frequency domain output using running DFT.
                        Data stored in HDF file output.h5 at location /floatingType/<name>
                        
     RCS_FDOM_HDF - bbox is surface or volume (in which case surfaces are selected according to mask).
                    Total power flow out of each surface of bbox determined.
                    Far field RCS is determined using far-field transform.
                    Frequency domain output using running DFT.
                    Data stored in HDF file output.h5 at location /floatingType/<name>

- Current/impedance observers?

- Add FDOM_BINARY arb. bbox/fields (derive from flim?)
                  used as utility by far-field transform!
                  write frequency.dat format?  

- Port observers.

  OP <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT COAXIAL     <r: innerRadius>/<r: Zc> <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]
  OP <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT CIRCULAR    <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]
  OP <ilo> <ihi> <jlo> <jhi> <klo> <khi> PORT RECTANGULAR <t: mode> [ <t: wfName> ] [ <r: size> ] [ <r:delay> ]

  Bounding box must be surface.
  Outer port size taken from bounding box.
  Modes: "TEM" "TE:1,0", "TM:1,1",....
  Use mode orthognality to determine mode "amplitude/phase" - voltage?
  Determine field propagation in one direction.
  Need to normalise relative to input port ampltide to get scattering parameters.

- Add NF-FF transformation.

  Do Poynting calculation on boundary and total raidated power 
  calculation on far-field.

  Geometric mean stuff?

### Add impedance observer for voltage sources (or generic)

Need to add support to observers to have internal observers that are 
used as services by FF, VR. internal attribute causes them to only
update internal state during main time loop. They don't dump but
their data is access later.

* TDOM_ASCII: find values and write to disk immediately.
* TDOM_BINARY: find values and write to disk immediately.
* FDOM_ASCII: find value and update internal array, dumped at end
* FDOM_IMPED: Create FDOM_ASCII obervers which behave except they 
  don't dump at end. Post-process data in observers.
* FDOM_FF: Create FDOM_ASCII obervers which behave except they don't 
  dump at end. Post-process data in observers.

If want proper averaging to nodes/times must cache at least last field values....

Geometric mean approach?

Use latest falcon implementation.

### Make observers report fields at bbox nodes

    #ifdef NODAL_OBSERVER
   
    inline real nodalEx( int i , int j , int k )
    {
      real zeta;
      zeta = dex[i-1] / ( dex[i] + dex[i-1] );
      return ( 1 - zeta ) * UNSCALE_Ex( Ex[i-1][j][k] , i - 1 ) + zeta * UNSCALE_Ex( Ex[i][j][k] , i );
    }

    inline real nodalEy( int i , int j , int k )
    {
      real zeta;
      zeta = dey[j-1] / ( dey[j] + dey[j-1] );
      return ( 1 - zeta ) * UNSCALE_Ey( Ey[i][j-1][k] , j - 1 ) + zeta * UNSCALE_Ey( Ey[i][j][k] , j );
    }

    inline real nodalEz( int i , int j , int k )
    {
      real zeta;
      zeta = dez[k-1] / ( dez[k] + dez[k-1] );
      return ( 1 - zeta ) * UNSCALE_Ez( Ez[i][j][k-1] , k - 1 ) + zeta * UNSCALE_Ez( Ez[i][j][k] , k );
    }

    inline real nodalHx( int i , int j , int k )
    {
      real zeta;
      real xi;
      zeta = dey[j-1] / ( dey[j] + dey[j-1] );
      xi = dez[k-1] / ( dez[k] + dez[k-1] );
      return ( 1 - zeta ) * ( 1 - xi ) * UNSCALE_Hx( Hx[i][j-1][k-1] , i ) +
             ( 1 - zeta ) * xi         * UNSCALE_Hx( Hx[i][j-1][k]   , i ) +
             zeta         * ( 1 - xi ) * UNSCALE_Hx( Hx[i][j][k-1]   , i ) +
             zeta         * xi         * UNSCALE_Hx( Hx[i][j][k]     , i );
    }

More complex - make observers report fields a evenly sampled physical positions 
within bounding box - use 3D interpolation to determine spatial averages.


## Internal boundaries

### Hacking SIBC to be compatible with PMC/PERIODIC BCs is very crude

Try and improve performance.

### Make SIBC PFE same as Debye

### Investigate late time instability of SIBC

### Check SIBC stability against changing Couraunt Number

### Check orientation and angle for isotropic symmetric case

### Test closed PEC box and closed SIBC box

## Check for ripple on time domain output from SIBC

### Small aperture subcell model

    BT <t: name> APERTURE <r: alpha_mxx> <r: alpha_myy>  <r: alpha_ez> <r: area> 
    TB <int: XLO> <int: XHI> <int: YLO> <int: YHI> <int: ZLO> <int: ZHI> \
       <t: name> [ <i: orient> [ <r: align> ] ]

Bounding box must be single face, e.g:

    BT circ_hole APERTURE 1e-6 1e-6 1e-6 1e-6

Implemented as effective surface impedance using first order model. Area needed 
to verify small aperture cut-off frequency. Electric polarisability not used,
requires more extensive changes. 


## External boundaries

### Check Mur interoperability with PML

If it doesn't work and too hard to fix disallow:

    if( thereAreExternalSurfaces( BT_MUR ) && thereAreExternalSurfaces( BT_PML ) )
      message( MSG_ERROR , 0 , "*** External PML and Mur boundary types cannot be used together.\n" );

### Add third order Higdon ABC with damping term

### Replace UPML with CFS-PML

Maybe not difficult since alot of the code with be the same. Make it compile time 
switchable `USE_CFS_PML=ON/OFF`.


## Wires

### Add thin wires

    WT <tag> THIN_WIRE radius
             PEC
    TW <xlo> <xhi> <ylo> <yhi> <zlo> <zhi> <tag> ....

Gedney notes show method using material coefficients only.

### Delta-gap voltage source

### End capacitance model

### Corners

### Resistive wires
