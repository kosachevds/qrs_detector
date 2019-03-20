#ifndef QRS_DETECTION_H
#define QRS_DETECTION

#define MARK_NO_QRS (0)
#define MARK_QRS    (1)

int DetectQrsPeaks(double const* signal, int size, char* result, double rate);

#endif