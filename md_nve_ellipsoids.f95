!##############################################################################!
! Programm zur Molekulardynamiksimulation von Ellipsoiden                      !
! periodische Randbedingungen, Würfel als Domain,                              !
! verwendet Perram-Wertheim-Kontaktfunktion                                    !
! liest ein file (inp.cnf) als Anfangskonfiguration ein. Format:               !
! Teilchenzahl                                                                 !
! Boxgröße                                                                     !
! r(3), Ausrichtung(3), v(3), omega(3)        (pro Teilchen eine Zeile)        !
! Gibt allgemeine Simulationsparameter aus                                     !
! kontrolliert auf Overlap in der Startkonfiguration                           !
!##############################################################################!

program md_nve_ellipsoids
  use var                  ! Variablen
  use ellipsoid_md         ! Sämtliche Funktionen und Subroutines
  implicit none

  ! Simulationsparameter
  real    :: box           ! Länge der Simulationsdomain (Würfel)
  real*8  :: time          ! Zeit nach welcher cnf geschrieben wird
  real*8  :: simtime

  ! Sonstige Variablen
  logical :: overlap       ! wenn .true. --> Überlapp
  integer :: i, j, k       ! Ellipsoidindex, Zähler für Stöße

  real :: start_time, stop_time     ! Zur Zeitnehmung
  real*8 :: starte, ende

  call cpu_time(start_time)         ! Speichere Beginnzeit

  open(1, file = "input")           ! Öffne Configurationsfile
  read(1, *) ha, hb, time, simtime  ! Lies erste Parameter ein
  

  ! Vorschläge für Simulationsparameter
  filecount = 1            ! nötig, um outputfiles richtig zu nennen
  ! time      = 1d0          ! Zeit, nach der cnf.xxx File geschrieben wird

  ! Ellipsenparamete
  ! ha  = 1.                 ! Hauptachsen in x und y Richtung (ungedreht)
  ! hb  = 2.                 ! Hauptachse in z Richtung (ungedreht)
  xab = ha/hb              ! Achsenverhältnis
  t_now =   0d0            ! Anfangszeit
  t_since = 0d0            ! Zeit seit letztem Stoß

  ! Simulationsparameter ausgeben
  write(*, '(a, t26, f15.4)') "Time between measurement", time

  ! Anfangskonfiguration einlesen
  call readcnf(box)

  ! auf Überlapp in Startkonfiguration prüfen
  call overlapcheck(overlap, i, j)
  if (overlap) then
    write(*, '(a, i3, i3, a)') 'Error: Overlap in startconfiguration', i, j
    stop
  else   
    write(*, '(a)') 'No overlap in startconfig - start simulation' 
    write(*, * )
  end if

  starte = calce(box)
  write(*, '(a, f25.10)') "Startenergy E: ", starte    ! Schreibe Startenergie


  do k = 1, 100
    ! call progress(1d0)
    call writecnf(k, box)
    call tocoll(i, j, box, time)
    call collide(box, i, j)
    print *,  "collision", k, i, j, t_now

    ! print *,  "T= ", t_now, " ", i, " ", j
    ! print *,  "E= ", calce(box)
    ! print *,  "----------------------------------------------------------------------"
  end do

  ! k = 0                                                ! Stoßzähler
  ! do while (t_now < simtime)                           ! Simuliere bis Simtime erreicht
  !   call tocoll(i, j, box, time)                       ! Lasse Zeit bis nächstem Überlapp voranschreiten
  !   do
  !     call collide(box, i, j)                          ! Berechne Stoß
  !     k = k+1                                          ! Zähle Stoß
  !     print *, "collision", k, i, j, t_now             ! Gib Stoß aus
  !     call overlapcheck(overlap, i, j, i, j+1)         ! Berechne ob noch weitere Teilchen stoßen
  !     if (overlap .eqv. .false.) exit                  ! Wenn keine weiteren Stöße wieder progress bis nächstem Stoß
  !   end do
  ! end do
  k = k-1
  write(*, '(a,i17.6)') "collisions calculated: ", k

  write(*, '(a, f29.8)') "Endtime T= ", t_now

  ende = calce(box)
  write(*, '(a, f27.10)') "Endenergy E: ", ende
  write(*, '(a, f29.10)') "Deltae dE: ", abs(ende-starte)

  call cpu_time(stop_time)
  print *,  
  write(*, '(a, f22.10)') "runtime (seconds):", stop_time - start_time

end program md_nve_ellipsoids

