module ellipsoid_md
  implicit none
  private
  public :: overlapcheck, tocoll, collide, readcnf, calce, progress, calcrotax, writecnf, mintc
  contains

!##############################################################################!
! Berechnet den Kontaktpunkt von zwei Ellipsoiden in box=1 Koordinaten
! berechnet und speichert auch den normierten Normalvektor
!##############################################################################!
  subroutine calccontact(i, j, rc, normal)
    use var, only : r, lambda
    implicit none
    integer                :: i, j            ! Teilchenindizes
    real*8, dimension(3)   :: rc              ! Ort des Stoßes
    real*8, dimension(3)   :: normal, rab     ! Normalvektor, Abstandsvektor
    real*8, dimension(3,3) :: mati, matj, mat ! Matrizen der Teilchen, kombinierte Matrix

    mati = ellmat(i)                          ! Berechne Teilchenmatrizen
    matj = ellmat(j)

    !mat = (1d0-lambda)*matj + lambda*mati     ! Berechne kombinierte Matrix
    mat = lambda*matj + (1d0-lambda)*mati

    mat = invert(mat)                         ! invertiere kombinierte Matrix

    rab(:) = r(:,j) - r(:,i) - anint(r(:,j) - r(:,i)) ! Abstandsvektor zw. i und j - PBC
    normal = matmul(mat, rab)                         ! Berechne Normalvektor

    !rc(:) = r(:,i) + lambda * matmul(mati,normal)   ! Berechne Kontaktpunkt
    rc(:) = r(:,i) + (1d0-lambda)*matmul(mati, normal)

    normal = normal/sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)  ! normiere Normalvektor
  end subroutine
       



  !##############################################################################!
  ! Berechnet die Normalgeschwindigkeit am Kontaktpunkt - siehe Stillinger
  !##############################################################################!
  function calcvn(i, j, rc, normal)
    use var, only : r, v, omega
    implicit none
    real*8                 :: calcvn          ! Normalgeschwindigkeit
    integer                :: i, j            ! Teilchenindizes
    real*8, dimension(3)   :: vcvec           ! vc als Vektor
    real*8, dimension(3)   :: rc, rac, rbc    ! Kontaktpunkt, Vektoren zum Kontaktpunkt
    real*8, dimension(3)   :: normal          ! normierter Normalvektor

    rac(:) = rc(:) - r(:,i) - anint(rc(:) - r(:,i))   ! Berechne Abstandsvektoren zu rc
    rbc(:) = rc(:) - r(:,j) - anint(rc(:) - r(:,j))   ! PBCs beachten

    ! Berechne vc
    vcvec(:) = (v(:,j) + cross(omega(:,j), rbc(:))) - &
            &  (v(:,i) + cross(omega(:,i), rac(:)))

    calcvn = dot_product(normal, vcvec)
    return
  end function calcvn


  !##############################################################################!
  ! Führe einen Stoß zwischen zwei Ellipsoiden durch (Stillinger)
  ! Wikipedia Collision Response
  !##############################################################################!
  subroutine collide(box, i, j)
    use var, only : r, omega, v, inertia, mass
    implicit none
    integer              :: i, j                     ! Teilchenindezes
    real*8               :: deltapab, vn             ! Impulsübertrag, Normalgeschwindigkeit
    real*8, dimension(3) :: normal, rc               ! normierter Normalvektor und Kontaktpunkt
    real*8, dimension(3) :: rca, rcb                 ! Abstandsvektoren zu rc
    real*8, dimension(3) :: cross1, cross2           ! Zwischenspeicher
    real :: box

    call calccontact(i, j, rc, normal)               ! Berechne Normalvektor und Kontaktpunkt
    vn = calcvn(i, j, rc, normal)                    ! Berechne Normalgeschwindigkeit

    if (vn .gt. 0) then
      print *,  "WARNING!!!!!!!!!!!!!!!!!!!!!! vn > 0!!!!!!!!!!!!!", i, j
      return
    end if

    !rca(:) = r(:,i) - rc(:) - anint(r(:,i) - rc(:))  ! Berechne Abstandsvektoren zu rc
    !rcb(:) = r(:,j) - rc(:) - anint(r(:,j) - rc(:))  ! PBCs beachten

    rca(:) = rc(:) - r(:,i) - anint(rc(:) - r(:,i))   ! eigentlich ist das rac
    rcb(:) = rc(:) - r(:,j) - anint(rc(:) - r(:,j))   ! eigentlich ist das rbc

    cross1 = cross(rca, normal)
    cross2 = cross(rcb, normal)

    !deltapab = 2d0*vn/ (1d0/mass + 1d0/mass + sqrt(dot_product(cross1, cross1))/ &    ! Berechne Impulsübertrag
    !                 & inertia + sqrt(dot_product(cross2, cross2))/inertia)

    vn = vn*box                                         ! Transformiere in sigma = 1 Koordinaten
    cross1 = cross1*box                                 ! ist nötig weil Impulsübertrag sonst falsch berechnet wird
    cross2 = cross2*box
    rca = rca*box
    rcb = rcb*box
    v = v*box

    deltapab = -2d0 * vn/ (1d0/mass + 1d0/mass + dot_product((1d0/inertia *  &   ! Impulsübertrag laut Wikipedia
            & cross(cross1, rca) + 1d0/inertia* &
            & cross(cross2, rcb)), normal))

    v(:,i) = v(:,i) - (deltapab/mass * normal(:))         ! Berechne neue Geschwindigkeit der Teilchen
    v(:,j) = v(:,j) + (deltapab/mass * normal(:))

    omega(:,i) = omega(:,i) - deltapab/inertia * cross1   ! Berechne neue Rotationsgeschwindigkeiten der Teilchen
    omega(:,j) = omega(:,j) + deltapab/inertia * cross2

    v = v/box                                             ! transformiere Geschwindigkeit zurück, alles andere wird verworfen
    vn = calcvn(i, j, rc, normal)                         ! unnötig?

    if (vn < 0) then
            
    print *,  vn
    end if
     
  end subroutine 



!##############################################################################!
! Funktion um die gesamte (kinetische) Energie zu berechnen
!##############################################################################!
  function calce(box)
    use var, only: omega, inertia, mass, v, n
    implicit none
    real*8  :: calce
    real    :: box
    integer :: i
    v = v*box
    calce = 0d0
    do i = 1, n
      calce = calce + (inertia * dot_product(omega(:,i), omega(:,i)) + &
                          mass * dot_product(v(:,i), v(:,i)))  * 5d-1
    end do

    v = v/box
    return
  end function calce


  function gettc(ell1, ell2, tmax)
    use var, only: r, axis
    implicit none

    integer :: ell1, ell2
    real*8, dimension(3,3) :: mat1, mat2      ! Aus Hauptachsen berechnete Matrizen

    real*8, dimension(3)     :: vr1, vr2        ! Ortsvektoren
    real*8, dimension(3)     :: vrab            ! Abstandsvektor zw. ell1 und ell2
    real*8                   :: fmin
    real*8, dimension(3)     :: startr1, startr2
    real*8, dimension(3,3)   :: startax1, startax2

    real*8 :: fab1
    real*8 :: fab2
    real*8 :: fabm
    real*8 :: t1
    real*8 :: t2
    real*8 :: tm
    real*8 :: lambda
    real*8 :: tmax
    real*8 :: gettc
    real*8 :: dt
    real*8 :: deltat
    integer :: k

    gettc = 0d0
    k = 0
    dt = 1d-4
    deltat = 1d0

    startr1 = r(:, ell1)
    startr2 = r(:, ell2)
    startax1 = axis(:,:,ell1)
    startax2 = axis(:,:,ell2)

    do
      k = k+1
      t2 = k*dt
      call movetwo(t2, ell1, ell2)
      mat1 = ellmat(ell1)
      mat2 = ellmat(ell2)
      vr1 = r(:, ell1)
      vr2 = r(:, ell2)
      vrab = vr2 - vr1 - anint(vr2 - vr1)
      call brent(0., 0.3, 1., 0.00001, fmin, mat1, mat2, vrab, lambda)
      fab2 = -pwcontact(lambda, mat1, mat2, vrab) - 1

      r(:, ell1) = startr1
      r(:, ell2) = startr2
      axis(:,:,ell1) = startax1
      axis(:,:,ell2) = startax2
      if (fab2 < 0) then 

        do while (deltat > 1d-10)
          tm = (t1+t2)/2
          call movetwo(tm, ell1, ell2)
          mat1 = ellmat(ell1)
          mat2 = ellmat(ell2)
          vr1 = r(:, ell1)
          vr2 = r(:, ell2)
          vrab = vr2 - vr1 - anint(vr2 - vr1)
          call brent(0., 0.3, 1., 0.00001, fmin, mat1, mat2, vrab, lambda)
          fabm = -pwcontact(lambda, mat1, mat2, vrab) - 1
          if (fab1*fabm < 0) then
            t2 = tm
          else
            t1 = tm
          end if
          r(:, ell1) = startr1
          r(:, ell2) = startr2
          axis(:,:,ell1) = startax1
          axis(:,:,ell2) = startax2

          deltat = abs(t1-t2)
        end do
        gettc = t1

        exit
      end if
      if (gettc > tmax) then
        gettc = 1d2
        exit
      end if
      gettc = gettc + dt
      fab1 = fab2
      t1 = t2
    end do
    r(:, ell1) = startr1
    r(:, ell2) = startr2
    axis(:,:,ell1) = startax1
    axis(:,:,ell2) = startax2
   
   return
  end function

  subroutine mintc(t, ell1, ell2)
    use var, only: n
    implicit none
    integer :: ell1, ell2
    real*8 :: t
    real*8 :: temp
    integer :: i, j

    t = 1d2
    do i = 1, n-1
      do j = i+1, n
        temp = gettc(i, j, t)
        if (temp < t) then
          t = temp
          ell1 = i
          ell2 = j
        end if
      end do
    end do
  end subroutine


  function calcnormal(mat1, mat2, vrab, lambda)
    implicit none
    real*8, dimension(3,3) :: mat1, mat2      ! Aus Hauptachsen berechnete Matrizen
    real*8, dimension(3,3) :: mat
    real*8, dimension(3)   :: vrab
    real*8, dimension(3)   :: calcnormal
    real*8                 :: lambda

    mat = lambda*mat2 + (1d0-lambda)*mat1

    mat = invert(mat)                         ! invertiere kombinierte Matrix
    calcnormal = matmul(mat, vrab)                         ! Berechne Normalvektor
    
  end function



!##############################################################################!
! Subroutine um auf Überlapp zu untersuchen                                    !
! berechnet aus den Achsenrichtungen für ein Teilchen (und alle nachkommenden) !
! die benötigten Matrizen um mit Perram/Wertheim Kontaktfunktion               !
! auf Überlapp zu überprüfen                                                   !
! übergibt an pwcontact die Matrizen um Kontaktfunktion zu berechnen           !
! Kontaktfunktion nur, wenn Abstand kleiner als doppelte Halbachsenlänge       !
!##############################################################################!
  subroutine overlapcheck(overlap, ell1, ell2, start, start2)
    use var, only: r, ha, hb, n, lambda
    implicit none

    logical :: overlap                        ! wenn true dann gibt es einen overlap
    integer :: i, j                           ! Schleifenvariablen
    integer :: ell1, ell2                     ! Indizes der Überlappenden Ellipsoide

    integer, optional, value :: start, start2 ! Indizes der ersten Ellipsoide die Überprüft werden
                                              ! werden diese Indezes nicht übergeben werden alle überprüft

    real*8, dimension(3,3) :: mati, matj      ! Aus Hauptachsen berechnete Matrizen

    real*8, dimension(3)   :: vec             ! Abstandsvektor zw. i und j
    real*8                 :: rab             ! Abstand zw. i und j
    real*8                 :: fmin

    real*8 :: one                             ! zum Vergleichen

    overlap = .false.                         ! Ausgegangen wird von keinem Überlapp
    one = 0.9999999d0

    if(.not. present(start)) then             ! Wenn Startindizes nicht übergeben werden 
      start  = 1                              ! überprüfe alle Ellipsoide
      start2 = 2
    end if

    if(start2 .eq. n+1) then                  ! Schutz vor falscher Parameterübergabe
      start  = start + 1
      start2 = start + 1
    end if

    ! gehe alle Teilchen ab start durch bis auf das letzte, weil immer mit Nachfolgern kontrolliert wird
    do i = start, n-1

      mati = ellmat(i)       ! Dyadisches Produkt für Matrix für Perram-Wertheim Kontaktfunktion
      
      do j = start2, n          ! Kontrolliere auf Überlapp von Teilchen i nur mit nachfolgenden Teilchen ab start2
        vec(:) = r(:,j) - r(:,i) - anint(r(:,j) - r(:,i))          ! Abstandsvektor zw. i und j - PBC
  
        rab = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))  ! Abstand zw. i und j
  
        if (rab .le. (max(ha,hb)*2)) then                          ! Überlapp nur möglich wenn rab kleiner als doppelte Halbachse

          matj = ellmat(j)   ! Dyadisches Produkt für Matrix für Perram-Wertheim Kontaktfunktion
  
    ! Jetzt hat man alles für die Kontaktfunktion. Bei gleichen Teilchen 
    ! (wie es hier der Fall ist, wenn die Ellipsoide ungedreht sind) 
    ! liegt das Maximum der Kontaktfunktion bei 0.5,
    ! für f == 1 --> Die Teilchen berühren sich
  
          ! ist der Wert der Kontaktfunktion an 0.5 größer als 1 liegt kein Überlapp vor
          ! wenn nicht, dann muss man noch weiter suchen, weil gedrehte Ellipsoide sind unterschiedliche Teilchen
          if (-pwcontact(5d-1, mati, matj, vec) .lt. one) then
                                                                       ! Suche noch weiter nach dem Minimum der Kontaktfunktion
            call brent(0., 0.3, 1., 0.00001, fmin, mati, matj, vec, lambda)
            fmin = -fmin

            if(lambda .lt. 1. .and. lambda .gt. 0.) then
              if(fmin .lt. one) then
                overlap = .true.
                ell1 = i            ! Indizes der Überlappenden Matrizen speichern
                ell2 = j

              end if
            end if

          end if

        end if
        if(overlap .eqv. .true.) exit  ! wenn schon Überlapp gefunden wurde muss nicht weiter kontrolliert werden
      end do
      start2 = i+2                     ! im ersten Durchlauf wird ab j=start2 überprüft, ab dann mit allen Nachfolgern
      if(overlap .eqv. .true.) exit
    end do
  end subroutine


!##############################################################################!
! Subroutine, die die Kontaktfunktion minimiert und das lambda(x), und den
! Funktionswert auf lambda und fmin speichert              
!##############################################################################!
  subroutine brent(ax, bx, cx, tol, fmin, MatI, MatJ, vec1, lambda)

  implicit none
  
    INTEGER :: ITMAX
  
    REAL :: ax,bx,cx,tol,CGOLD,ZEPS
    real*8 :: fmin
    PARAMETER (ITMAX=50,CGOLD=.3819660, ZEPS=1.0e-10)
    real*8,dimension (3,3) :: MatI,MatJ
    real*8,dimension (3) :: vec1 
    real*8 :: lambda
  
  !Given a function f, and given a bracketing triplet of abscissas ax,bx,cx
  !(such that bx is between ax and cx and f(bx) is greater than both f(ax) and
  !f(cx) ), this  routine  isolates the maximum to a fractional precision of about
  !tol using Brent's method.  The abscissa of the  maximum  is  returned as
  !xmin, and  the maximum  function  value  is  returned as brent.
  
  !Parameters:  Maximum allowed number of iterations; golden ratio
  
    INTEGER:: iter
    REAL*8:: a,b,d,e,etemp,fu,fv,fw,fx,p,q,k,tol1,tol2,u,vv,w,x,xm
  
    a = min(ax,cx)
    b = max(ax,cx)
    vv = bx
    w = vv
    x = vv
    e = 0.                          !This will be the distance moved on the step before last.
    fx = pwcontact(x, MatI, MatJ, vec1)
    fv = fx
    fw = fx
  
  
     do iter = 1,ITMAX
   
        xm = 0.5*(a+b)
        tol1 = tol*abs(x)+ZEPS
        tol2 = 2.*tol1
  
        if(abs(x-xm) .le. (tol2-.5*(b-a))) goto 3     !Test for done  here.
  
        if(abs(e) .gt. tol1) then                     !Construct a trial  parabolic fit
  
           k = (x-w)*(fx-fv)
           q = (x-vv)*(fx-fw)
           p = (x-vv)*q-(x-w)*k
           q = 2.*(q-k)
  
           if (q .gt. 0.) p = -p
  
           q = abs(q)
     etemp = e
     e = d
  
     if (abs(p) .ge. abs(.5*q*etemp) .or. p .le. q*(a-x) .or. p .ge. q*(b-x)) goto 1
  
  !The above conditions determine the acceptability of the parabolic fit
  !Here it is o.k.:
  
     d = p/q                    !Take  the  parabolic  step.
     u = x+d
  
     if (u-a .lt. tol2 .or. b-u .lt. tol2) d = sign( tol1, xm-x )
  
     goto 2                    !Skip  over the  golden section step.
  
        endif
   
  1     if (x .ge. xm) then                    !We  arrive  here  for a  golden  section  step, which  we  take
                        !into  the larger  of the two segments.
           e = a-x
  
        else
        
           e = b-x
  
        endif
  
        d = CGOLD*e                     !Take  the  golden  section step.
  
  2     if (abs(d) .ge. tol1) then                    !Arrive here with d computed either from parabolic fit, or
                               !else  from  golden section.
     u = x+d
  
        else
  
           u = x+sign(tol1,d)
  
        endif
  
        fu = pwcontact( u, MatI, MatJ, vec1 )                    !This is  the one function evaluation per iteration,
  
        if (fu .lt. -1.) then                    !die kontaktfunktion bei überlapp bleibt < 1, lieg kein überlapp vor
                      !ist das minimum hier < -1 --> fertig
           x = u
     fx = fu
     goto 3
  
        end if
  
        if (fu .le. fx) then                               !and now we have to decide what to do with our function
                                       !evaluation.  Housekeeping follows:
     if (u .ge. x) then
      
        a = x
    
     else
  
        b = x
  
     endif
  
     vv = w
     fv = fw
     w = x
     fw = fx
     x = u
     fx = fu
  
        else
  
     if(u .lt. x) then
  
        a = u
  
     else
  
        b = u
  
     endif
   
     if(fu .le. fw .or. w .eq. x) then
  
        vv = w
        fv = fw
        w = u
        fw = fu
  
     else if(fu .le. fv .or. vv .eq. x .or. vv .eq. w) then
  
        vv = u
        fv = fu
  
     endif
  
        endif
  
  !Done with housekeeping.  Back for another iteration.
     enddo
  
  !   pause 'brent exceed maximum iterations'
  
  3   lambda = x  !Arrive  here ready  to exit  with  best values.
      fmin = fx
  
  end subroutine





!##############################################################################!
! Funktion um Matrix für eine Ellipse zu berechnen (X_A - siehe Stillinger)    !
!##############################################################################!
  function ellmat(i)
    use var, only : axis
    integer                :: k, l     ! Schleifenvariable
    integer                :: i        ! Teilchenindex
    real*8, dimension(3,3) :: ellmat   ! Matrix der Ellipse
    ! Dyadisches produkt für Matrix für Perram-Wertheim Kontaktfunktion
    do k = 1, 3
      do l = 1, 3
        ellmat(l,k) = axis(1,l,i)*axis(1,k,i) + axis(2,l,i)*axis(2,k,i) + &
                &   axis(3,l,i)*axis(3,k,i)
      end do
    end do
    return
  end function ellmat




!##############################################################################!
! Subroutine um für Teilchen i die Hauptachsen zu berechnen.                   !
! (achsennummer, Richtung)                                                     !
!##############################################################################!
  subroutine ha_calc(i)
    use var, only : axis, ha, hb!, angle
    implicit none
    integer                :: i             ! Teilchenindex von dem Hauptachsen berechnet werden
    !real*8, dimension(3)   :: cosi, sini    ! Zwischenspeicher für Winkelfunktionen

    !sini(1) = sin(angle(1,i))               ! Werden öfter gebraucht darum gleich gespeichert
    !sini(2) = sin(angle(2,i))
    !sini(3) = sin(angle(3,i))

    !cosi(1) = cos(angle(1,i))
    !cosi(2) = cos(angle(2,i))
    !cosi(3) = cos(angle(3,i))

    !! Siehe Wikipedia - Euler Rotation zxz Konvention
    !axis(1,1,i) = cosi(1)*cosi(3) - sini(1)*cosi(2)*sini(3)     ! Achse a
    !axis(1,2,i) = sini(1)*cosi(3) + cosi(1)*cosi(2)*sini(3)
    !axis(1,3,i) = sini(2)*sini(3)
    !!axis(1,1,i) =  cosi(1)*cosi(3) - sini(1)*cosi(2)*sini(3)   ! Achse a original?? war da die Drehung falsch???
    !!axis(1,2,i) = -cosi(1)*sini(3) - sini(1)*cosi(2)*cosi(3)
    !!axis(1,3,i) =  sini(1)*sini(2)
  
    !axis(2,1,i) = -cosi(1)*sini(3) - sini(1)*cosi(2)*cosi(3)
    !axis(2,2,i) = -sini(1)*sini(3) + cosi(1)*cosi(2)*cosi(3)
    !axis(2,3,i) =  sini(2)*cosi(3)
    !!axis(2,1) =  sini(1)*cosi(3) + cosi(1)*cosi(2)*sini(3)   ! Achse a
    !!axis(2,2) = -sini(1)*sini(3) + cosi(1)*cosi(2)*cosi(3)
    !!axis(2,3) =  cosi(1)*sini(2)
  
    !axis(3,1,i) =  sini(1)*sini(2)
    !axis(3,2,i) = -cosi(1)*sini(2)
    !axis(3,3,i) =  cosi(2)
    !!axis(3,1,i) = sini(2)*sini(3)                              ! Achse b
    !!axis(3,2,i) = sini(2)*cosi(3)
    !!axis(3,3,i) = cosi(2)

    !axis(1,:,i) = axis(1,:,i) * ha          ! Achsen auf Richtige Länge bringen
    !axis(2,:,i) = axis(2,:,i) * ha
    !axis(3,:,i) = axis(3,:,i) * hb
    axis(2,:,i) = cross(axis(3,:,i),axis(1,:,i))
    axis(2,:,i) = axis(2,:,i)/norm(axis(2,:,i))*ha
    axis(1,:,i) = axis(1,:,i)/norm(axis(1,:,i))*ha
    axis(3,:,i) = axis(3,:,i)/norm(axis(3,:,i))*hb
  end subroutine 



!##############################################################################!
! Funktion die aus den Ellipsenmatrizen, die in overlapcheck berechnet werden  !
! den Wert der Kontaktfunktion berechnet
!##############################################################################!
  real*8 function pwcontact(x, mati, matj, vec)
    implicit none
  
    real*8, dimension(3,3) :: mat, matinv     ! Matrix + Inverse aus Ellipsenmatrizen
    real*8, dimension(3,3) :: mati, matj      ! Die Matrizen der Ellipsen
    real*8, dimension(3)   :: vec, vec2       ! Speicher für Zwischenergebnisse
    real*8                 :: x               ! lambda

    
    mat = x*mati + (1d0-x)*matj  

    matinv = invert(mat)       ! invertiere Matrix

    vec2 = matmul(matinv,vec)  ! Kontaktfunktion berechnen
    pwcontact = -x*(1d0-x) * dot_product(vec,vec2)

    return
  end function pwcontact


!##############################################################################!
! Funktion um Inverse einer 3x3 Matrix zu berechnen
!##############################################################################!
  function invert(mat)
    real*8, dimension(3,3) :: mat, invert
    real*8, parameter      :: zero = 1d-10    ! Null, um Regularität zu prüfen
    real*8                 :: det             ! Determinante

    ! Berechne Determinante
    det = mat(1,1)*mat(2,2)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1) + &
        & mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1) - &
        & mat(2,3)*mat(3,2)*mat(1,1) - mat(3,3)*mat(1,2)*mat(2,1)

    ! Prüfe ob singuläre Matrix vorliegt
    if (abs(det) .le. zero) then
      write(*, '(a)') 'Singular matrix'
      stop
      return
    end if

    ! berechne inverse
    invert(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
    invert(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
    invert(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
    invert(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
    invert(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
    invert(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
    invert(3,1) = mat(2,1)*mat(3,2) - mat(3,1)*mat(2,2)
    invert(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
    invert(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

    invert = invert / det

    return
  end function invert




!##############################################################################!
! Subroutine um die Teilchen in der angegebenen Zeit t zu bewegen              !
! (Translation + Rotation), periodische RB werden berücksichtigt               !
!##############################################################################!
  subroutine progress(t)
    use var, only : r, omega, v, n
    implicit none
    real*8               :: t               ! So viel Zeit soll verstreichen
    integer              :: i 
    real*8               :: temp_angle      ! Drehwinkel, aus Winkelgeschwindigkeit berechnet


    r(:,:)     = r(:,:) + t * v(:,:)        ! Translation (box = 1 Einheiten)
    r(:,:)     = r(:,:) - anint(r(:,:))     ! Periodische Randbedingungen

    ! t_now = t_now + t

    ! Rotation der Ellipsoide
    do i = 1, n
      temp_angle = norm(omega(:,i)) * t     ! Berechne Drehwinkel aus Winkelgeschwindigkeit
      if (temp_angle .eq. 0d0) cycle
      call quatrot(i, temp_angle)           ! Rotiere jedes Teilchen mittels Quaternionen
    end do
  end subroutine 

  subroutine movetwo(t, i, j)
    use var, only : r, omega, v
    implicit none
    real*8  :: t
    integer :: i, j
    real*8  :: temp_angle

    r(:,i)     = r(:,i) + t * v(:,i)        ! Translation (box = 1 Einheiten)
    r(:,j)     = r(:,j) + t * v(:,j)        ! Translation (box = 1 Einheiten)
    r(:,i)     = r(:,i) - anint(r(:,i))     ! Periodische Randbedingungen
    r(:,j)     = r(:,j) - anint(r(:,j))     ! Periodische Randbedingungen

    temp_angle = norm(omega(:,i)) * t     ! Berechne Drehwinkel aus Winkelgeschwindigkeit
    if (temp_angle .ne. 0d0) call quatrot(i, temp_angle)
    temp_angle = norm(omega(:,j)) * t
    if (temp_angle .ne. 0d0) call quatrot(i, temp_angle)
  end subroutine

!##############################################################################!
! Funktion die die Norm eines 3D Vektors berechnet
!##############################################################################!
  function norm(vect)
    implicit none
    real*8, dimension(3) :: vect
    real*8               :: norm
    norm = sqrt(dot_product(vect, vect))
    return
  end function norm


!############################################################################!
! Dreht Teilchen i um rotangle um seine drei Achsen                          !
! verwendet dabei Quaternionen                                               !
!############################################################################!
  subroutine quatrot(i, rotangle)
    use var, only : axis, omega!, ha, hb
    implicit none
    integer                :: i, j                ! Teilchenindex, Schleifenvariable
    real*8                 :: rotangle            ! Drehwinkel um jede Achse
    !real*8, dimension(3,4) :: quater              ! enthält für jede Drehachse ein Quaternion
    real*8, dimension(4)   :: tempquat            ! Zwischenspeicher
    real*8, dimension(4)   :: resquat, iresquat   ! resultierendes Quaternion aus Kombination von den drei Drehungen
    real*8, dimension(3)   :: direction           ! normalisierter Richtungsvektor in Drehachsenrichtung
    real*8, dimension(4)   :: axisquat            ! Quaternion in Richtung der Achse
    !real*8, dimension(4)   :: omegaquat
  
    ! Berechne die drei Quaternionen für die einzelnen Drehungen
    !do j = 1, 3
    !  direction(:) = axis(j,:,i) / &
    !          & sqrt(axis(j,1,i)**2 + axis(j,2,i)**2 + axis(j,3,i)**2)
    !  quater(j,:) = calcquat(rotangle(j), direction(:))
    !end do


    ! Kombiniere die einzelnen Drehungen zu einer kombinierten
    !tempquat =  quatproduct(quater(1, :), quater(2, :))
    !resquat  =  quatproduct(tempquat  , quater(3, :))

    resquat(1) = cos(rotangle/2.)
    direction(:) = omega(:,i) / norm(omega(:,i))
    resquat(2:4) = direction(:) * sin(rotangle/2)

    iresquat(1)   = resquat(1)                  ! Invertiertes Rotationsquaternion
    iresquat(2:4) = -resquat(2:4)

    ! Rotiere die drei Achsen
    do j = 1, 3
      axisquat(1)   = 0d0                         ! Setze Achsenquaternion
      axisquat(2:4) = axis(j, :, i)
      tempquat = quatproduct(iresquat, axisquat)  ! Rotiere die Achse
      axisquat = quatproduct(tempquat, resquat)

      axis(j, :, i) = axisquat(2:4)               ! Lies die Achse aus dem Quaternion
    end do

    !axis(1,:,i) = axis(1,:,i) / norm(axis(1,:,i)) * ha  ! Achsen renormieren
    !axis(2,:,i) = axis(2,:,i) / norm(axis(2,:,i)) * ha
    !axis(3,:,i) = axis(3,:,i) / norm(axis(3,:,i)) * hb
  end subroutine 



!############################################################################!
! Funktion die für ein Teilchen i, die Drehachse berechnet
! und in Omega (Betrag & Richtung) speichert
!############################################################################!
  subroutine calcrotax(i)
    use var, only : axis, omega
    implicit none
    real*8, dimension(3) :: rotax
    integer :: i, j
    real*8, dimension(3)   :: direction           ! normalisierter Richtungsvektor in Drehachsenrichtung
    real*8, dimension(3,4) :: quater              ! enthält für jede Drehachse ein Quaternion
    real*8, dimension(4)   :: tempquat            ! Zwischenspeicher
    real*8, dimension(4)   :: resquat             ! resultierendes Quaternion aus Kombination von den drei Drehungen

    do j = 1, 3
      direction(:) = axis(j,:,i) / &
              & sqrt(axis(j,1,i)**2 + axis(j,2,i)**2 + axis(j,3,i)**2)
      quater(j,:) = calcquat(omega(j,i), direction(:))
    end do

    ! Kombiniere die einzelnen Drehungen zu einer kombinierten
    tempquat =  quatproduct(quater(1, :), quater(2, :))
    resquat  =  quatproduct(tempquat  , quater(3, :))

    if (sin(acos(resquat(1))) .eq. 0) then
        rotax(:) = 0d0
        omega(:,i) = 0d0
        return
    end if

    rotax(:) = resquat(2:4)/ sin(acos(resquat(1)))
    rotax(:) = rotax(:) / sqrt(rotax(1)**2 + rotax(2)**2 + rotax(3)**2)

    omega(:,i) = sqrt(dot_product(omega(:,i), omega(:,i))) * rotax(:)
  end subroutine




!############################################################################!
! Berechnet ein Quaternion für einen Bestimmten Winkel und Drehachse         !
!############################################################################!
  function calcquat(angle, direction)
    implicit none
    real*8               :: angle       ! Drehwinkel
    real*8, dimension(3) :: direction   ! Achsenrichtung
    real*8, dimension(4) :: calcquat    ! Quaternion
    integer              :: i           ! Schleifenvariable

    calcquat(1) = cos(angle/2)
    do i = 2, 4
      calcquat(i) = sin(angle/2) * direction(i-1)
    end do
    return
  end function calcquat



  !############################################################################!
  ! Berechnet das Quaternionenprodukt von zwei Quaternionen                    !
  !############################################################################!
  function quatproduct(quat1, quat2)
    implicit none
    real*8, dimension(4) :: quat1, quat2, quatproduct
    real*8, dimension(3) :: c
  
    c =  cross(quat1(2:4), quat2(2:4))       ! Kreuzprodukt vom Vektoranteil

    quatproduct(1)   = quat1(1)*quat2(1)   - dot_product(quat1(2:4),quat2(2:4)) 
    quatproduct(2:4) = quat1(1)*quat2(2:4) + quat2(1)*quat1(2:4) - c
    return
  end function quatproduct



  !############################################################################!
  ! Berechnet das Kreuzprodukt von zwei 3D Vektoren                            !
  !############################################################################!
  function cross(a, b)
    implicit none
    real*8, dimension(3) :: a, b
    real*8, dimension(3) :: cross

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
    return
  end function cross




  !############################################################################!
  ! Bewegt die Teilchen bis zur nächsten Kollision                             !
  ! Schreibt in regelmäßigen Abständen die momentane Konfiguration             !
  !############################################################################!
  subroutine tocoll(i, j, box, dt)
    use var, only : r, axis, filecount, n, t_now, t_since
    implicit none
    real*8, dimension(3,n)   :: startr       ! Startorte der Teilchen
    real*8, dimension(3,3,n) :: startaxis    ! Startwerte der Axen
    logical                  :: overlap      ! wenn .true. -> Überlapp
    integer                  :: i, j         ! Indizes der Stoßenden Teilchen
    integer                  :: k, l, m      ! Zähler der Zeitschritte
    real*8                   :: bigt, smallt ! grober und feiner Zeitschritt
    real*8                   :: xst          ! sehr kleiner Zeitschritt
    real*8                   :: coltime      ! Zeit bis zum nächsten Stoß
    real                     :: box          ! Simulationsbox
    real*8                   :: dt           ! Zeitschritt, nach welchem cnf ausgegeben werden soll
    ! ----------------------------------
    integer :: a, b
    ! ----------------------------------

    ! Initialisierungen
    bigt    = 1d-2
    smallt  = 1d-4
    xst     = -1d-5
    overlap = .false.
    k = 0
    l = 0
    m = 0
    ! speichere Anfangszustand
    startr(:,:)      = r(:,:)
    startaxis(:,:,:) = axis(:,:,:)

    ! berechne grob die Zeit bis Zusammenstoß
    do while (overlap .eqv. .false.)
      call progress(bigt)
      call overlapcheck(overlap, i, j)
      k = k+1
      if (k .gt. 10000) exit
    end do

    ! Setze Teilchen zurück
    r(:,:)      = startr(:,:)
    axis(:,:,:) = startaxis(:,:,:) 

    ! Wenn nach weniger als 1000 groben Zeitschritten ein Überlapp zustandekam -> verfeinere Zeitmessung
    if (overlap .eqv. .true.) then
      coltime = bigt * (k - 1)               ! Zeit bis knapp vor Überlapp

      overlap = .false.
      call progress(coltime)                 ! unechter Progress bis knapp vor Überlapp

      do while (overlap .eqv. .false.)       ! unechter Progress bis Stoß
        call progress(smallt)
        call overlapcheck(overlap, i, j)
        l = l+1
      end do



      r(:,:)      = startr(:,:)              ! Setze Teilchen zurück an den Start
      axis(:,:,:) = startaxis(:,:,:) 

      coltime = bigt*(k-1) + smallt*l        ! Verfeinerte Kollisionszeit

      ! -----------------------------------------------
      call progress(coltime)
      do while (overlap .eqv. .true.)
        a = i
        b = j
        call progress(xst)
        m = m + 1
        call overlapcheck(overlap, i, j)
      end do
      i = a
      j = b
      coltime = coltime + m*xst
      ! -----------------------------------------------

      r(:,:)      = startr(:,:)              ! Setze Teilchen zurück an den Start
      axis(:,:,:) = startaxis(:,:,:) 

      do while (coltime + t_since .ge. dt)   ! echter Progress  mit dazwischen Fileausgabe
        call progress(dt-t_since)            ! bis nächstem schreiben
        coltime = coltime - (dt-t_since)     ! Ziehe von Kollisionszeit die schon vergangene Zeit ab
        t_now = t_now + (dt-t_since)
        call writecnf(filecount, box)
        filecount = filecount + 1
        t_since = 0d0
      end do

      call progress(coltime)             ! restliche Zeit bis zur Kollision
      t_since = t_since + coltime
      t_now = t_now + coltime

    else                                 ! Wenn kein Stoß passiert
      print *,  l, "l"
      print *,  k, "k"
      write(*, '(a)') 'Error: No next collision'
      stop
    end if
  end subroutine



  !############################################################################!
  ! Subroutine um Anfangskonfiguration einzulesen                              !
  ! allokiert auch Speicher                                                    !
  ! transformiert von sigma = 1 Koordinaten in Boxlänge = 1 Koordinaten        !
  ! berechnet Anfangsorientierung der Hauptachsen                              !
  !############################################################################!
  subroutine readcnf(box)
    use var, only : r, v, ha, hb, omega, axis, n, t_now
    implicit none
    integer           :: i, fileerror                     ! Schleifenvariable, Filecheck
    real              :: box                              ! Domaingröße
    character (len=400) :: tmp
  
    ! Teilchenzahl und Domaingröße einlesen und ausgeben
    open(100, file = "inp.cnf")
    read(100, fmt = *, iostat = fileerror) n
    if (fileerror /= 0) then
        write(*, *) "Error while reading inp.cnf"
        stop
    end if
    read(100, fmt = *, iostat = fileerror) box
    if (fileerror /= 0) then
        write(*, *) "Error while reading inp.cnf"
        stop
    end if
    read(100, *) t_now
    read(100, *) tmp
    write(*, '(a, t26, i15)')    "Number of Particles", n
    write(*, '(a, t26, f15.4)')  "Boxsize",             box
    write(*, '(a, t26, f15.4)')  "ha",                  ha
    write(*, '(a, t26, f15.4)')  "hb",                  hb

    ! Allocate arrays 
    allocate(r(3,n), v(3,n), omega(3,n), axis(3,3,n))

    ! Teilchenpositionen, Ausrichtungen und Geschwindigkeiten einlesen
    do i = 1, n
      read(100, *, iostat = fileerror) r(:,i), axis(3,:,i), axis(1,:,i), v(:,i), omega(:,i)
      if (fileerror /= 0) then
        write(*, *) "Error while reading inp.cnf"
        stop
      end if
    end do

    r(:,:) = r(:,:) / box                     ! Transformiere auf box = 1
    r(:,:) = r(:,:) - anint(r(:,:))           ! Periodische RB
    v(:,:) = v(:,:) / box
    ha = ha / box                             ! Transformiere auf box = 1
    hb = hb / box

    do i = 1, n                              
      call ha_calc(i)                         ! Berechne die Hauptachsen für alle Teilchen
      call calcrotax(i)                       ! Berechne richtiges Omega für alle Teilchen
    end do


    close(100)
    write(*, '(a)') "read startconfiguration"
    call writecnf(0, box)
  end subroutine 





  !############################################################################!
  ! Subroutine um Konfiguration zu schreiben                                   !
  ! Format:                                                                    !
  ! t_now                                                                      !
  ! r(3), za(3), xa(3), v(3), omega(3)        (pro Teilchen eine Zeile)        !
  ! Transformiert auch gleich von box = 1 Koordinaten in Sigma = 1 Koordinaten !
  !############################################################################!
  subroutine writecnf(outcnf, box)
    use var, only : r, v, t_now, omega, n, axis
    implicit none
    integer          :: outcnf, i         ! Dateinummer, Schleifenvariable
    character(len=7) :: filename          ! Ausgabefiles: cnf.xxx
    real             :: box               ! Simulationsdomainlänge

    r(:,:) = r(:,:) * box                 ! Transformiere auf sigma = 1 Koordinaten
    v(:,:) = v(:,:) * box
    axis(:,:,:) = axis(:,:,:) * box

    ! Öffne File mit richtigem Namen
    write(filename, "('cnf.', I3.3)") outcnf
    open(10, file = filename)

    write(10, "(i4.4)"),  n
    write(10, "(f6.2)"),  box
    write(10, "(f10.5)"), t_now ! Schreibe aktuelle Zeit
    write(10, "(a2, a13, a13, a12, a12, a12, a12, a12, a12, a11, a12, a12, a16, a12, a12)") &
            & "rx", "ry", "rz", "zax", "zay", "zaz", "xax", "xay", "xaz", "vx", "vy", "vz", "omegax", "omegay", &
            & "omegaz"
  
    ! Schreibe Teilchendaten
    do i = 1, n
      write(10, "(f11.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,&
         & f12.8,f12.8,f12.8,f12.8,f12.8,f12.8,f12.8)") r(:,i), axis(3,:,i), axis(1,:,i), v(:,i), omega(:,i)
    end do
    close(10)

    r(:,:) = r(:,:) / box                 ! Transformiere zurück auf box = 1
    v(:,:) = v(:,:) / box
    axis(:,:,:) = axis(:,:,:) / box
  end subroutine 

end module ellipsoid_md
