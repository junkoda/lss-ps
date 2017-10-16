

for dix in range(3):
    if dix == 0:
        print('if(ix_left <= ix0[0] && ix0[0] < ix_right) {')
    elif dix == 1:
        print('if(ix_left <= ix1[0] && ix1[0] < ix_right) {')
    elif dix == 2:
        print('if(ix_left <= ix2[0] && ix2[0] < ix_right) {')
        
    for diy in range(3):
        for diz in range(3):
            print('    d->add(ix%d[0], ix%d[1], ix%d[2], w*w%d[0]*w%d[1]*w%d[2]);' % (dix, diy, diz, dix, diy, diz))

    print('}')
