
program recode_genotypes

implicit none

integer :: nsnp,io,j,anim,nlines,ns
integer(kind=1), allocatable, dimension(:) :: SNP
integer(kind=4), allocatable, dimension(:) :: ind
integer(kind=1), allocatable, dimension(:,:) :: allele
character(len=35) :: format_5,format_6
character(len=200) :: fline

!Read in index that stores segregating SNPs
open(1,file='segregating_loci.txt')
ns=0
do
 read(1,*,iostat=io)j
 if(io/=0)exit
 ns=ns+1
enddo
rewind(1)
allocate(ind(ns))
do j=1,ns
 read(1,*)ind(j)
enddo
close(1)

!Count number of SNPs ("nsnp")
nsnp = -1       !Account for first line
open(3,file="lm_mrk_001.txt")
do
  read(3,*,iostat=io)
  if ( io /= 0 ) exit
  nsnp = nsnp + 1
end do
close(3)
print *,nsnp,' SNPs were detected in file lm_mrk_001.txt'

open(5,file="p1_mrk_001.txt")
!Format: ID      Genotypes (0 = a1,a1; 2 = a2,a2; 3 = a1,a2; 4 = a2,a1; 5 = missing; The first allele is paternal and the second allele is maternal) ...
read(5,*)       !skip the first line
read(5,'(a)')fline
call derive_format_statement(nsnp,format_5,fline)
call derive_format_statement(ns,format_6,fline)
backspace(5)
open(6,file="phased_alleles.txt")
open(7,file="SNP_genotypes.txt")
allocate(SNP(nsnp))
allocate(allele(nsnp,2))
nlines = 0
do
  read(5,format_5,iostat = io)anim,SNP
  if ( io /= 0 ) exit
  !recode to get two lines of phased alleles
  do j = 1, nsnp
   if ( SNP(j) == 0 ) then
    allele(j,1:2) = 0 !a1 = 0
   else if ( SNP(j) == 2 ) then
    allele(j,1:2) = 1 !a2 = 1
   else if ( SNP(j) == 3 ) then
    allele(j,1) = 0
    allele(j,2) = 1
   else if ( SNP(j) == 4 ) then
    allele(j,1) = 1
    allele(j,2) = 0
   end if  
  end do
  write(6,format_6)anim,allele(ind(1:ns),1)
  write(6,format_6)anim,allele(ind(1:ns),2)
  !now recoded to get genotypes on 0,1,2 scale
  do j = 1, nsnp
    if ( SNP(j) > 2 ) SNP(j) = 1
  end do
  write(7,format_5)anim,SNP(ind(1:ns))
  nlines = nlines + 1
end do
close(5)
close(6)
close(7)
print *,nlines,' animals were detected in file p1_mrk_001.txt.'

contains
subroutine derive_format_statement(n_snp,fformat,fline)
implicit none
integer :: n_snp,rep_30,max_nr,rest_nr,n_snp2
integer, dimension(3) :: nr_spaces      !stores: 1) number of preceding spaces; number of digits for the ID; number of spaces between ID and genotypes
 character(len=*) :: fformat
 character(len=*) :: fline

call count_spaces(fline,nr_spaces)

max_nr = 30000

if ( nr_spaces(1) == 0 .and. nr_spaces(3) == 0 ) then
        if ( n_snp > max_nr ) then
          rep_30 = int(n_snp/max_nr)
          fformat = '(i07,0000(30000i1),00000i1)'
          rest_nr = n_snp - rep_30 * max_nr
          write(fformat(6:9),'(i4.4)')rep_30
          write(fformat(20:24),'(i5.5)')rest_nr
        else
        !  fformat = '(00x,i00,00x,00000i1)'
          fformat = '(i07,00000i1)'
          write(fformat(6:10),'(i5.5)')n_snp
        end if
else
        if ( n_snp > max_nr ) then
          rep_30 = int(n_snp/max_nr)
        !  fformat = '(00x,i00,00x,00(30000i1),00000i1)'
          fformat = '(i00,00x,0000(30000i1),00000i1)'
          rest_nr = n_snp - rep_30 * max_nr
          write(fformat(10:13),'(i4.4)')rep_30
          write(fformat(24:28),'(i5.5)')rest_nr
        else
        !  fformat = '(00x,i00,00x,00000i1)'
          fformat = '(i00,00x,00000i1)'
          write(fformat(10:14),'(i5.5)')n_snp
        end if
        !write(fformat(2:3),'(i2.2)')nr_spaces(1)
        write(fformat(3:4),'(i2.2)')sum(nr_spaces(1:2))
        write(fformat(6:7),'(i2.2)')nr_spaces(3)
end if
!print *,' The genotype file is read using format statement ',fformat

end subroutine derive_format_statement

subroutine count_spaces(sentence,nr_spaces)
!counts number of words in sentence
!a word is defined as a substring consisting of one or more characters that are not blanks

INTEGER :: nr_spaces(3), kk
 CHARACTER(LEN=*), INTENT(IN) :: sentence

INTEGER :: I, stringlength
INTEGER,ALLOCATABLE,DIMENSION(:) :: space
 CHARACTER :: charac

kk = nr_words(sentence)
! if ( kk /= 2 ) then
  ! print *,' Problem reading genotypes without spaces between columns.'
  ! print *,' The file should contain one separate column with IDs, followed by genotypes.'
  ! print *,sentence
  ! stop
! end if

stringlength = LEN(sentence)

allocate(space(0:stringlength)); space = 0
do i=1,stringlength
   charac = sentence(I:I)
   if ( charac .EQ. ' ' ) then
     if ( space(i-1) < 2 ) then !a space was detected on the previous position, or this position wasn't evaluated yet (space(i-1)=0)
       space(i) = 1                             !a space was detected on position i, and this is before the ID.
         else if ( space(i-1) <= 3 ) then
       space(i) = 3
         else
       space(i) = 9                             !this makes sure that possible spaces after the genotypes (within the first 100 positions), are added to the spaces between ID & genotypes.
         end if
   else
     if ( i == 1 ) then                 !the ID starts on the first position (i.e. there is no space on the first position).
       space(i) = 2
         else
       if ( space(i-1) <= 2 ) then
         space(i) = 2
           else
             space(i) = 9
           end if
         end if
   end if
!   print *,i,space(i)
end do

!determine number of spaces (1), length of ID (2), and number of spaces between ID and genotypes (3)
do i = 1, 3
  nr_spaces(i) = count(space == i)
end do

if ( nr_spaces(1) > 1 ) then
  nr_spaces(2) = nr_spaces(2) + nr_spaces(1) - 1
  nr_spaces(1) = 1
end if

deallocate(space)

end subroutine count_spaces

!*...............................................................................
function nr_words(sentence)
!counts number of words in sentence
!a word is defined as a substring consisting of one or more characters that are not blanks

INTEGER:: nr_words
CHARACTER(LEN=*), INTENT(IN) :: sentence

INTEGER :: I, stringlength
CHARACTER :: lastchar

stringlength = LEN(sentence)
nr_words = 0
lastchar = ' '

do i=1,stringlength
   if ( (lastchar .EQ. ' ') .AND. (sentence(I:I) .NE. ' ') ) then
      ! a blank followed by another character indicates the start of a new word
      nr_words = nr_words + 1
   end if
   lastchar = sentence(I:I)
end do

end function nr_words
!*......................................................................................

end program recode_genotypes
