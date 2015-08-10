import numpy as np
import cv2

cap = cv2.VideoCapture("D:\\data\\wmv_files\\140314-XY5-H2B dendra- hRSV-RFP.wmv")

while(cap.isOpened()):
    ret, frame = cap.read()

    #gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    gray = frame[:,:,2]
    cv2.imshow('frame',gray)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

cap.release()
cv2.destroyAllWindows()