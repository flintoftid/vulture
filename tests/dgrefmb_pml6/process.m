
opNum = 1;
timeNum = 1:270;
freqNum = 1:27;

cd( '../dgref_ref' );

[ xB , yB , zB , exB , eyB , ezB , hxB , hyB , hzB , tB ] = tdfdReadImpulseDat3D( opNum , timeNum );
[ xB , yB , zB , ExB , EyB , EzB , HxB , HyB , HzB , fB ] = tdfdReadFrequencyDat3D( opNum , freqNum , 1 );

cd( '../dgref_pml6' );

[ xT , yT , zT , exT , eyT , ezT , hxT , hyT , hzT , tT ] = tdfdReadImpulseDat3D( opNum , timeNum );
[ xT , yT , zT , ExT , EyT , EzT , HxT , HyT , HzT , fT ] = tdfdReadFrequencyDat3D( opNum , freqNum , 1 );

LRE_ex = exT - exB;
LRE_ey = eyT - eyB;
LRE_ez = ezT - ezB;

LRE_Ex = ExT - ExB;
LRE_Ey = EyT - EyB;
LRE_Ez = EzT - EzB;

maxLRE_ex = max( abs( LRE_ex(:) ) );
maxLRE_ey = max( abs( LRE_ey(:) ) );
maxLRE_ez = max( abs( LRE_ez(:) ) );

maxLRE_Ex = max( abs( LRE_Ex(:) ) );
maxLRE_Ey = max( abs( LRE_Ey(:) ) );
maxLRE_Ez = max( abs( LRE_Ez(:) ) );

LRE_e2 = abs( LRE_ex ).^2 + abs( LRE_ey ).^2 + abs( LRE_ez ).^2;

LRE_E2 = abs( LRE_Ex ).^2 + abs( LRE_Ey ).^2 + abs( LRE_Ez ).^2;

GRE_e = squeeze( sum( sum( sum( LRE_e2 , 1 ) , 2 ) , 3 ) );  

GRE_E = squeeze( sum( sum( sum( LRE_E2 , 1 ) , 2 ) , 3 ) ); 

maxGRE_e = max( GRE_e(:) );

maxGRE_E = max( GRE_E(:) );

CGRE_e = sum( GRE_e );

CGRE_E = sum( GRE_E );

fprintf( 'Max. LRE ex = %.1f dB\n' , db20( maxLRE_ex ) );
fprintf( 'Max. LRE ey = %.1f dB\n' , db20( maxLRE_ey ) );
fprintf( 'Max. LRE ez = %.1f dB\n' , db20( maxLRE_ez ) );
fprintf( 'Max. GRE e  = %.1f dB\n' , db10( maxGRE_e ) );
fprintf( 'Cum, GRE e  = %.1f dB\n' , db10( CGRE_e ) );

fprintf( 'Max. LRE Ex = %.1f dB\n' , db20( maxLRE_Ex ) );
fprintf( 'Max. LRE Ey = %.1f dB\n' , db20( maxLRE_Ey ) );
fprintf( 'Max. LRE Ez = %.1f dB\n' , db20( maxLRE_Ez ) );
fprintf( 'Max. GRE E  = %.1f dB\n' , db10( maxGRE_E ) );
fprintf( 'Cum, GRE E  = %.1f dB\n' , db10( CGRE_E ) );

figure( 1 );
hl1 = plot( tT / 1e-9 , db10( GRE_e ) , 'r-o' );
hxl = xlabel( 'Time (ns)' );
hyl = ylabel( 'Global reflection error (dB)' );
print( 'dgrefmb_pml6_GRE_e' , '-depsc2' );

figure( 2 );
hl1 = plot( fT / 1e9 , db10( GRE_E ) , 'r-o' );
hxl = xlabel( 'Frequency (GHz)' );
hyl = ylabel( 'Global reflection error (dB)' );
print( 'dgrefmb_pml6_GRE_E' , '-depsc2' );

figure( 3 );
fidx = 10;
fGHz =  fT(10) / 1e9;
hl1  = plot( xT , db10( squeeze( LRE_E2(:,1,1,fidx) ) ) , 'k-o;YLO-ZLO;' );
hold on;
hl2  = plot( xT , db10( squeeze( LRE_E2(:,1,7,fidx) ) ) , 'r-o;YLO-ZHI;' );
hl3  = plot( xT , db10( squeeze( LRE_E2(:,7,1,fidx) ) ) , 'g-o;YHI-ZLO;' );
hl4  = plot( xT , db10( squeeze( LRE_E2(:,7,7,fidx) ) ) , 'b-o;YHI-ZHI;' );
hl5  = plot( yT , db10( squeeze( LRE_E2(1,:,1,fidx) ) ) , 'k-*;XLO-ZLO;' );
hl6  = plot( yT , db10( squeeze( LRE_E2(1,:,7,fidx) ) ) , 'r-*;XLO-ZHI;' );
hl7  = plot( yT , db10( squeeze( LRE_E2(7,:,1,fidx) ) ) , 'g-*;XHI-ZLO;' );
hl8  = plot( yT , db10( squeeze( LRE_E2(7,:,7,fidx) ) ) , 'b-*;XHI-ZHI;' );
hl9  = plot( zT , db10( squeeze( LRE_E2(1,1,:,fidx) ) ) , 'k-^;XLO-YLO;' );
hl10 = plot( zT , db10( squeeze( LRE_E2(1,7,:,fidx) ) ) , 'r-^;XLO-YHI;' );
hl11 = plot( zT , db10( squeeze( LRE_E2(7,1,:,fidx) ) ) , 'g-^;XHI-YLO;' );
hl12 = plot( zT , db10( squeeze( LRE_E2(7,7,:,fidx) ) ) , 'b-^;XHI-YHI;' );
hxl = xlabel( 'i, j ,k (cells)' );
hyl = ylabel( 'Local reflection error (dB)' );
hti = title( sprintf( 'Frequency %.1f GHz' , fGHz ) );
print( sprintf( 'dgrefmb_pml6_LRE_E2_%.1fGHz.eps' , fGHz ) , '-depsc2' );

assert( db20( maxLRE_ex ) < -125 );
assert( db20( maxLRE_ey ) < -125 );
assert( db20( maxLRE_ez ) < -125 );
assert( db10( maxGRE_e )  < -105 );
assert( db20( maxLRE_Ex ) < -38 );
assert( db20( maxLRE_Ey ) < -38 );
assert( db20( maxLRE_Ez ) < -38);
assert( db10( maxGRE_E )  < -15 );

