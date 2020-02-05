# Copyright: Luis Pedro Coelho <luis@luispedro.org>, 2012
# License: MIT
import numpy as np

def read_roi(fileobj):
    '''
    points = read_roi(fileobj)

    Read ImageJ's ROI format
    '''
# This is based on:
# https://urldefense.proofpoint.com/v2/url?u=http-3A__rsbweb.nih.gov_ij_developer_source_ij_io_RoiDecoder.java.html&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=dhePJ2vCFP-Wcc6a2dOYN_oh0vldYxUfyaCE0cDBqf4&m=xjDE8rCSVszUIsfq_iJR50c_qeN3jje5P32lE5LW5oQ&s=qonm2SFzcPKawr7mtdBZBEDtAubm_jxwA2pa8Fj-YJw&e= 
# https://urldefense.proofpoint.com/v2/url?u=http-3A__rsbweb.nih.gov_ij_developer_source_ij_io_RoiEncoder.java.html&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=dhePJ2vCFP-Wcc6a2dOYN_oh0vldYxUfyaCE0cDBqf4&m=xjDE8rCSVszUIsfq_iJR50c_qeN3jje5P32lE5LW5oQ&s=2EDrL0LCIff2TFxfTA7Zr-usfr6uec3enhlJlTWKDQY&e= 


    SPLINE_FIT = 1
    DOUBLE_HEADED = 2
    OUTLINE = 4
    OVERLAY_LABELS = 8
    OVERLAY_NAMES = 16
    OVERLAY_BACKGROUNDS = 32
    OVERLAY_BOLD = 64
    SUB_PIXEL_RESOLUTION = 128
    DRAW_OFFSET = 256


    pos = [4]
    def get8():
        pos[0] += 1
        s = fileobj.read(1)
        if not s:
            raise IOError('readroi: Unexpected EOF')
        return ord(s)

    def get16():
        b0 = get8()
        b1 = get8()
        return (b0 << 8) | b1

    def get32():
        s0 = get16()
        s1 = get16()
        return (s0 << 16) | s1

    def getfloat():
        v = np.int32(get32())
        return v.view(np.float32)

    magic = fileobj.read(4)
    if magic != 'Iout':
        raise IOError('Magic number not found')
    version = get16()

    # It seems that the roi type field occupies 2 Bytes, but only one is used
    roi_type = get8()
    # Discard second Byte:
    get8()

    if not (0 <= roi_type < 11):
        raise ValueError('roireader: ROI type %s not supported' % roi_type)

    if roi_type != 7 and roi_type != 3 and roi_type != 5:
        raise ValueError('roireader: ROI type %s not supported (!= 7,3)' % roi_type)

    top = get16()
    left = get16()
    bottom = get16()
    right = get16()
    n_coordinates = get16()

    x1 = getfloat() 
    y1 = getfloat() 
    x2 = getfloat() 
    y2 = getfloat()
    stroke_width = get16()
    shape_roi_size = get32()
    stroke_color = get32()
    fill_color = get32()
    subtype = get16()
    if subtype != 0:
        raise ValueError('roireader: ROI subtype %s not supported (!= 0)' % subtype)
    options = get16()
    arrow_style = get8()
    arrow_head_size = get8()
    rect_arc_size = get16()
    position = get32()
    header2offset = get32()
    
    if(roi_type == 3):
        points = np.array([x1,y1,x2,y2])
      
        #roi.setDrawOffset(drawOffset); ??     
      
    elif(roi_type == 7 or roi_type == 5):

        if options & SUB_PIXEL_RESOLUTION:
            getc = getfloat
            points = np.empty((n_coordinates, 2), dtype=np.float32)
        else:
            getc = get16
            points = np.empty((n_coordinates, 2), dtype=np.int16)
        points[:,0] = [getc() for i in xrange(n_coordinates)]
        points[:,1] = [getc() for i in xrange(n_coordinates)]
        points[:,0] += left
        points[:,1] += top
        points -= 1
    return points

def read_roi_zip(fname):
    import zipfile
    with zipfile.ZipFile(fname) as zf:
        return [read_roi(zf.open(n))
                    for n in zf.namelist()]