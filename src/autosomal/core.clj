(ns autosomal.core
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]))

(defn read-gzipped-csv [input]
  (with-open [in (io/reader
                  (java.util.zip.GZIPInputStream.
                   (io/input-stream
                    input)))]
    (->> (csv/read-csv in)
         rest
         (map #(assoc {} :ch (get % 1) :pos (Long/parseLong (get % 2))))
         doall)))

(defn write-gzipped-csv [input filename]
  (with-open [out (io/writer
                   (java.util.zip.GZIPOutputStream.
                    (io/output-stream filename)))]
    (csv/write-csv out input)))

(def clean-dna
  (->> (read-gzipped-csv "randi.csv.gz")
       (remove #(= "0" (:ch %)))
       (remove #(= "X" (:ch %)))))

(defn filter-last [ch]
  (:pos (last (filter #(= ch (:ch %)) clean-dna))))

(def chromo-size37 (into [] (map filter-last (map str (range 1 23)))))

(defn random-dna-letter []
  (case (rand-int 4)
    0 "A"
    1 "C"
    2 "G"
    3 "T"))

(defn random-dna []
  (reduce #(conj %1 (assoc %2 :gen (random-dna-letter))) [] clean-dna))

(def father {1 (random-dna) 2 (random-dna)})

(def mother {1 (random-dna) 2 (random-dna)})

(defn randomize-crossovers-2 [snps-count]
  (let [num-of-crossovers (inc (rand-int 4))]
    (conj
     (into
      []
      (sort
       (for [_ (range num-of-crossovers)] (rand-int snps-count))))
     (dec snps-count))))

(defn randomize-crossovers [ch]
  (let [snps (into [] (sort-by :pos (filter #(= ch (:ch %)) clean-dna))) 
        rand-xover (randomize-crossovers-2 (count snps))]
    (map #(:pos (get snps %)) rand-xover)))

(defn generate-dna-from-parent [parent]
  (loop [parent-1 (get parent 1)
         parent-2 (get parent 2)
         result []
         ch "1"
         cx (randomize-crossovers ch)
         side (rand-int 2)]
    (if (seq parent-1)
      (let [p1 (first parent-1)
            p2 (first parent-2)]
        (if (= ch (:ch p1))
          (if (< (:pos p1) (first cx))
            (let [p (if (= 0 (mod side 2)) p1 p2)]
              (recur (rest parent-1) (rest parent-2) (conj result p) ch cx side))
            (let [p (if (= 0 (mod side 2)) p2 p1)]
              (recur (rest parent-1) (rest parent-2) (conj result p) ch (rest cx) (inc side))))
          (let [new-ch (:ch p1)
                new-cx (randomize-crossovers new-ch)
                new-side (rand-int 2)
                p (if (= 0 (mod side 2)) p1 p2)]
            (recur (rest parent-1) (rest parent-2) (conj result p) new-ch new-cx new-side))))
      result)))

(defn generate-child []
  (loop [dad (generate-dna-from-parent father)
         mum (generate-dna-from-parent mother)
         result []]
    (if (seq dad)
      (let [d (first dad)
            m (first mum)]
        (recur (rest dad) (rest mum) (conj result (assoc d :gen (str (:gen d) (:gen m))))))
      result)))

(defn make-vector [data]
  (reduce #(conj %1 [(:pos %2) (:ch %2) (:pos %2) (:gen %2)]) [["RSID" "CHROMOSOME" "POSITION" "RESULT"]] data))

(defn save-child [name]
  (-> (generate-child)
      (make-vector)
      (write-gzipped-csv name)))

(comment

  (sort-by :pos (filter #(= "2" (:ch %)) clean-dna))

  (save-child "test3.csv.gz")

  (generate-child)

  (zipmap [:a :b :c :d] (cycle [1 2]))

  (randomize-crossovers "1")

  (generate-dna-from-parent father)

  (filter #(= "0" (:ch %)) (get father 1))

  (filter #(= "X" (:ch %)) clean-dna)
  mother
  (random-dna)

  chromo-size37

  (last (filter #(= "1" (:ch %)) clean-dna))


  (write-gzipped-csv  [[1 2 3] [4 5 6] [7 8 9]] "test.csv.gz")
  )