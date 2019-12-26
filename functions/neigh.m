function [mov_col,mov_row]=neigh(dir)

switch dir
    case 1
        mov_col=1;
        mov_row=0;
    case 8
        mov_col=1;
        mov_row=-1;
    case 7
        mov_col=0;
        mov_row=-1;
    case 6
        mov_col=-1;
        mov_row=-1;
    case 5
        mov_col=-1;
        mov_row=0;
    case 4
        mov_col=-1;
        mov_row=+1;
    case 3
        mov_col=0;
        mov_row=+1;
    case 2
        mov_col=1;
        mov_row=+1;
end