
for FILE in animations/*.mp4; do
   NEW_FILE=${FILE// /_}
   echo ${NEW_FILE}
   echo ${FILE} ${NEW_FILE}
   mv "${FILE}" "${NEW_FILE}"
done


for FILE in animations/*.mp4; do 
   NAME=${FILE##*/}; 
   echo ${NAME}
   NEW_FILE=gifs/${NAME%.*}.gif
   echo ${NEW_FILE}
   echo ${FILE} ${NEW_FILE}
   ffmpeg -i ${FILE} ${NEW_FILE}
done
