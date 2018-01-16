
!##############################################################################!
! enthält die wichtigsten Variablen und arrays                                 !
!##############################################################################!
module var
  implicit none
  private
  
  ! public data
  integer,                                public :: n       ! Teilchenzahl
  real*8,  dimension(:,:),   allocatable, public :: r       ! Positionen (3,n)
  real*8,  dimension(:,:),   allocatable, public :: v       ! Geschwindigkeiten (3,n)
  real*8,  dimension(:,:),   allocatable, public :: angle   ! Ausrichtung (3,n) Eulerwinkel zxz-Konvention
  real*8,  dimension(:,:),   allocatable, public :: omega   ! Winkelgeschwindigkeit (3,n) (ha, ha, hb)
  real*8,  dimension(:,:,:), allocatable, public :: axis    ! Halbachsen (3,3,n) (Achse, Richtung, Teilchen)

  real, public :: xab       ! Achsenverhältnis
  real, public :: ha, hb    ! Ellipsenhauptachsen

  real, public, parameter :: mass    = 1.   ! Masse der Teilchen
  real, public, parameter :: inertia = 1.   ! Trägheitsmoment

  real*8,  public :: t_now     ! aktuelle Zeit
  real*8,  public :: t_since   ! Zeit seit letzter Ausgabe
  integer, public :: filecount ! Zähler für Ausgabefile

  real*8, public :: lambda     ! xwert des Maximums der Kontaktfunktion
end module var

