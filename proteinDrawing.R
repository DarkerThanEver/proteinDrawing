library(R6)
library(ggplot2)
library(dplyr)
library(drawProteins)
library(stringr)

# Two functions to semi-automatically create the three different tables needed
# Note: for all three it essential that the number of elements for the columns
#       as specified in the arguments is the same (as always for a data.frame)
# Note: if no arguments are specified in the function call then an empty data.frame
#       will be returned
createTheoryTable <- function(type = as.character(),
                              description = as.character(),
                              begin = as.numeric(),
                              end = as.numeric(),
                              length = as.numeric(),
                              accession = as.character(),
                              entryName = as.character(),
                              taxid = as.numeric(),
                              order = as.numeric()){
  return(data.frame(type = type,
                    description = description,
                    begin = begin,
                    end = end,
                    length = length,
                    accession = accession,
                    entryName = entryName,
                    taxid = taxid,
                    order = order,
                    stringsAsFactors = FALSE)
  )
}

createExperimentTable <- function(type = as.character(),
                                  description = as.character(),
                                  begin = as.numeric(),
                                  end = as.numeric(),
                                  length = as.numeric(),
                                  order = as.numeric(),
                                  sequence = as.character(),
                                  mods = as.character(),
                                  expMr = as.numeric(),
                                  calcMr = as.numeric(),
                                  rt = as.numeric(),
                                  expect = as.numeric(),
                                  varModPos = as.character()){
  return(data.frame(type = type,
                    description = description,
                    begin = begin,
                    end = end,
                    length = length,
                    order = order,
                    sequence = sequence,
                    mods = mods,
                    expMr = expMr,
                    calcMr = calcMr,
                    rt = rt,
                    expect = expect,
                    varModPos = varModPos,
                    stringsAsFactors = FALSE)
  )
}

chainData <- R6Class("chainData",
                     public = list(
                       table = NA,
                       # Essentially virtual function
                       # note: used to have descendents create correct table columns
                       createTable = function(){
                         self$Table <- data.frame()
                       },
                       # initialize data (deleteOrder = NA) or select data to wipe deleteOrder = 1 or 1:3 or c(1,3,4)
                       # self$Table MUST always contain a column called order 
                       wipeData = function(deleteOrder = NA){
                         if (identical(deleteOrder,NA)){
                           self$createTable()
                         } else {
                           self$table <- self$table %>% filter(!(order %in% deleteOrder))
                         }
                       },
                       # returns a data.frame with the data of only a certain "Type", eg "MOD_RES"
                       # select which proteins via orders, can also be done via the usual
                       # commands via dplyr etc
                       type = function(types = "MOD_RES", orders = NA){
                         orders <- unique(self$table$order)
                         return(self$table %>% filter(type %in% types & order %in% orders))
                       },
                       # takes a table (eg proteinFormatted mascot export table) and returns a data.frame with all
                       # the begin- and end-aa's; needs type = "COVER". Note: the end-aa is the first aa NOT in the
                       # specific protein section
                       # coverString: this can be used in the case of more than one sequence of coverage data in the table
                       #              eg when coverage data is present of different digestions
                       # note: by definition the function uses the self$table of the object, it can be
                       # note: duplicates with the same "sequence" are removed
                       findCoverage = function(data = self$table, coverString = "COVER", whichOrder = self$orders){
                         data <- data %>%
                           filter((type == coverString) & (order %in% whichOrder))
                         if ("sequence" %in% colnames(data)){
                           data <- data %>%
                             distinct(sequence, order, .keep_all = TRUE)  
                         }
                         return(data %>%
                                  dplyr::select(begin, end, order))
                       },
                       # finds all the modifications (modification) in a theoryTable or experimentTable
                       # like (TRUE/FALSE) allows for not using the exact name of the modification, it will
                       # eg allow Phospho and recognize both Phospho (ST) and Phospho (Y). Note however
                       # that this works via grepl, so using phospho (in this example) will NOT work
                       # outsideCoverage = TRUE serves to allow also mods that are not inside the coverage
                       #                   (=identified) parts of the protein. If set to FALSE it will check for this,
                       #                   otherwise not
                       # coverageTable: if outsideCoverage == FALSE and the table (features) is external (eg uniprot)
                       #                data, then include coverageTable (experimentTable) which contains type == "COVER"
                       #                and a column "sequence". Note: NOT the result from findCoverage, but an 'normal'
                       #                (experiment) table with rows of type == "COVER"
                       # coverString: this can be used in the case of more than one sequence of coverage data in the table
                       #              eg when coverage data is present of different digestions
                       findMods = function(data = self$table, type = "MOD_RES", modification, like = TRUE,
                                            outsideCoverage = TRUE,
                                            coverageTable = NA, coverString = "COVER"){
                         theData <- data
                         the_features_combined <- theData[-c(1:nrow(theData)),]
                         allOrders  <- unique(theData$order)
                         for (counter in 1:length(allOrders)){
                           data <- theData[(theData$type == type) & (theData$order == allOrders[counter]), ]
                           thelist <- as.logical()
                           if (!like){
                             for (counter2 in 1:nrow(data)){
                               thelist <- append(thelist,(data$description[counter2] == modification))
                             }
                           } else {
                             thelist <- grepl(modification,data$description)
                           }
                           the_features <- data[thelist, ]
                           if (nrow(the_features)>0){
                             if (!outsideCoverage){
                               doNotRemove <- rep(FALSE,nrow(the_features))
                               if (identical(coverageTable,NA)){
                                 takeCoverage <- self$findCoverage(theData[theData$order == allOrders[counter],],
                                                                   coverString = coverString)
                               } else {
                                 takeCoverage <- self$findCoverage(coverageTable[coverageTable$order == allOrders[counter],],
                                                                   coverString = coverString)
                               }
                               for (counter in 1:nrow(the_features)){
                                 for (counter2 in 1:nrow(takeCoverage)){
                                   if (the_features$begin[counter] %in% (takeCoverage$begin[counter2]:takeCoverage$end[counter2]) |
                                       the_features$end[counter] %in% (takeCoverage$begin[counter2]:takeCoverage$end[counter2])
                                   ){
                                     doNotRemove[counter] <- TRUE
                                   }
                                 }
                               }
                               the_features_combined <- bind_rows(the_features_combined,the_features[doNotRemove,])
                             } else {
                               the_features_combined <- bind_rows(the_features_combined,the_features)
                             }
                           }
                         }
                         return(the_features_combined)
                       },
                       # wrapper for draw_canvas, non essential
                       draw_canvas = function(data = self$table, ...){
                         return(drawProteins::draw_canvas(data, ...))
                       },
                       # function to add a graphics object onto an existing graph (p), data = eg,
                       # outInMin & outInMax specify how the positioning of the coverage bar, 
                       # coverString allows the use of multiple coverage data, eg different enzymes
                       # note: dependent on the "type" column containing rows with "COVER" in them
                       # note: use color = ... fill = ...  to specify the looks
                       draw_coverage = function (p, data = self$table, inOutMin = -0.20, inOutMax = 0.20, coverString = "COVER",...){
                         begin = end = description = NULL
                         data <- data %>%
                           filter(type == coverString)
                         p <- p + ggplot2::geom_rect(data = data[data$type == coverString,],
                                                     mapping = ggplot2::aes(xmin = begin, xmax = end-1,
                                                                            ymin = order + inOutMin, ymax = order + inOutMax),
                                                     ...)
                         return(p)
                       },
                       
                       # function to add a graphics object onto an existing graph (p), data = mascot theoryTable or experimentTable,
                       # outIn specifies the height where the modification should be place, uses findMods to search if the mod
                       # is present and where. See findMods for info on the parameters
                       # note: use color = ... fill = ...  size = ... shape = ... to specify the looks of the modifications
                       # note: a modification is defined by both the type and description columns. The mod argument = the description
                       draw_mods = function (p, data = self$table, modification = NA, type = "MOD_RES",
                                             like = TRUE, outsideCoverage = TRUE,
                                              coverageTable = NA, inOut = 0.25,...){
                         if (!identical(modification,NA)){
                           begin = end = description = NULL
                           p <- p + ggplot2::geom_point(data = self$findMods(data = data, type = type,
                                                                             modification = modification,
                                                                             like = like,
                                                                             outsideCoverage = outsideCoverage,
                                                                             coverageTable = coverageTable), 
                                                        ggplot2::aes(x = begin, y = order + inOut),...)
                         }
                         return(p)
                       },
                       # adds a line to the graphic object p based on the table
                       # the line is the connection between cysteine aa's when there is disulfide bridge
                       # dependent of "DISULFID" in the type column, and begin and end columns
                       # note: inOut determines the position of both the begining and the end of the line
                       # note: use color = ... size = ... , type = ...  to specify the looks of the line
                       draw_bridge = function(p, data = self$table, inOut = 0.25,...){
                         ordersInThere <- unique(data$order)
                         for (counter in 1:length(ordersInThere)){
                           data2 <- data[data$type == "DISULFID" & data$order == ordersInThere[counter],]
                           if (!identical(data2$type[1],as.character(NA))){
                             p <- p + geom_segment(data = data2,
                                                   aes(x=begin, xend=end,
                                                       y = order + inOut, yend = order + inOut), show.legend = FALSE, ...)
                           }
                         }
                         return(p)
                       },
                       # adds colored bars to the graphics object p depending on the 3D structure (helix strand or turn) as specified
                       # in the table data. The table will usually come from info from eg Swiss-/UniProt. This function is nearly identical
                       # to the original function in the package/library drawProteins; the only real change are the inOutMin & inOutMax
                       # parameters which allows specifying how much the colored bars are allowed to "stick out". A possible future change
                       # might be to add the ability to change the colors which are at the moment rather basic
                       draw_folding = function (p, data = self$table, show.legend = TRUE,
                                                 show_strand = TRUE, show_helix = TRUE, show_turn = TRUE,
                                                 inOutMin = -0.2, inOutMax = 0.2, ...){
                         begin = end = description = type = NULL
                         if (show_strand == TRUE) {
                           p <- p + ggplot2::geom_rect(data = dplyr::filter(data, 
                                                                            grepl("STRAND", type)), mapping = ggplot2::aes(xmin = begin, 
                                                                                                                           xmax = end, ymin = order + inOutMin, ymax = order + inOutMax, 
                                                                                                                           fill = type), show.legend = show.legend, ...)
                         }
                         if (show_helix == TRUE) {
                           p <- p + ggplot2::geom_rect(data = data[data$type == 
                                                                     "HELIX", ], mapping = ggplot2::aes(xmin = begin, 
                                                                                                        xmax = end, ymin = order + inOutMin, ymax = order + inOutMax,
                                                                                                        fill = type), show.legend = show.legend, ...)
                         }
                         if (show_turn == TRUE) {
                           p <- p + ggplot2::geom_rect(data = data[data$type == 
                                                                     "TURN", ], mapping = ggplot2::aes(xmin = begin, 
                                                                                                       xmax = end, ymin = order + inOutMin, ymax = order + inOutMax,
                                                                                                       fill = type), show.legend = show.legend, ...)
                         }
                         return(p)
                       },
                       # adds colored bars to the graphics object p depending on the data specified by type & description in the table data.
                       # originally meant to be able to show specific amino acid sequences in a chain
                       # if needed you can optionally do the selection of rows manually and leave both arguments empty
                       # modified function from drawProteins. inOutMin & inOutMax are used to position the graph objects (in-/outside
                       # the protein chain), use color = ... fill = ... to change the look of the bars
                       draw_seqPart = function (p, data = self$table, type =NA, description = NA, inOutMin = -0.2, inOutMax = 0.2,...) 
                       {
                         if (identical(type,NA)){
                           if (!identical(description,NA)){
                             data <- data[data$description == description,]
                           } # else do nothing, use all data, assumes manual data selection
                         } else {
                           if (identical(description,NA)){
                             data <- data[data$type == type,]
                           } else {
                             data <- data[(data$type == type) & (data$description == description),]
                           }
                         }
                         p <- p + ggplot2::geom_rect(data = data,
                                                     mapping = ggplot2::aes(xmin = begin,
                                                                            xmax = end,
                                                                            ymin = order + inOutMin,
                                                                            ymax = order + inOutMax),
                                                     ...)
                         return(p)
                       },
                       initialize = function(){
                         self$createTable()
                       }
                     ),
                     active = list(
                       orders = function(value){
                         if (missing(value)) {
                           return(unique(self$table$order))
                         } else {
                           # do nothing, read only
                         }
                       }
                     )
)

theoryChain <- R6Class("theoryChain",
                       inherit = chainData,
                       public = list(
                         createTable = function(){
                           self$table <- createTheoryTable()
                         },
                         # Function to pre-process the data.frame coming from either a feature_to_dataframe() call or
                         # a file containing the info in csv-format. Example: to get the data for ovalbumin (Gallus Gallus) 
                         # get_features("P01012") %>% feature_to_dataframe() %>% addInfo()
                         # theData: if dataMmode == FALSE then theData = a character string, a fileName that needs to be loaded
                         #          (the function assumes a csv formatted file)
                         # dataMode: if TRUE, then data comes from data.frame via get_features() %>% feature_to_dataframe(),
                         #           example: data for ovalbumin (Gallus Gallus)
                         #           get_features("P01012") %>% feature_to_dataframe() %>% readInfo()
                         #          if FALSE, then data comes from file, see theData
                         # whichInfo: specify which types of data are to be in the data.frame, if you leave whichInfo NA,
                         #            it will return all of the info, otherwise just the specified 'types'
                         #            Example: whichInfo = c("CARBOHYD","MOD_RES","DISULFID").
                         # correct: when using your own sequences (in eg Mascot searches), the actual numbers of the begin & end positions
                         #          may need updating. Example: when using Ovalbumin (Gallus Gallus) the first amino acid (Methionine) is
                         #          not present in the "finished" protein, so position 2 is actually the first amino acid!
                         #          --> use correct = 1
                         # order: in case you need to change the order of the protein in the data.frame for drawing
                         #        Note: this function will not work properly with data.frames containing more than one order
                         #        It will set all order values to the specified one
                         # rowNames: TRUE if dataMode == FALSE, then the first column will be removed (usually rownames), otherwise not
                         addTable = function(theData = NA, dataMode = TRUE, whichInfo = NA ,correct = 0, order = NA, rowNames = TRUE){
                           if (!identical(theData, NULL)){
                             if (!dataMode){
                               info_sp <- read.csv(theData, stringsAsFactors = FALSE, header = TRUE)
                               if (rowNames){
                                 info_sp <- info_sp[,-1]
                               }
                             } else {
                               info_sp <- theData
                             }
                             if (!identical(whichInfo,NA)){
                               info_sp <- info_sp %>% filter(type %in% whichInfo)
                             } 
                             info_sp[,3:4] <- info_sp[,3:4]-correct
                             if (!is.na(order)){
                               info_sp$order <- order
                             }
                             rownames(info_sp) <- NULL
                             self$table <- bind_rows(self$table, info_sp)
                           }
                         },
                         # same as addTheory, but now can work with either theData = list of filenames to load (dataMode = FALSE)
                         # or a data.frame containing more than one value for order
                         # note: correct = list of numbers to use for corrects of the different files/data
                         #       order = list of numbers to use for the order of the different files/data
                         #       whichInfo = as in addTheory, same for all different files/orders in the data.frame
                         #       dataMode = as in addTheory, same for all different files/orders in the data.frame
                         addTables = function(theData = NA, dataMode = TRUE, whichInfo = NA ,correct = NA, order = NA, rowNames = TRUE){
                           if (!identical(theData, NULL)){
                             if (!dataMode){
                               if (identical(correct,NA)){
                                 correct <- rep(0,length(theData))
                               }
                               if (identical(order,NA)){
                                 order <- rep(NA, length(theData))
                               }
                               for (counter in 1:length(theData)){
                                 self$addTable(theData[counter], dataMode = FALSE, whichInfo = whichInfo, correct = correct[counter],
                                                order = order[counter], rowNames = rowNames)
                               }
                             } else {
                               ordersInThere <- unique(theData$order)
                               if (identical(correct,NA)){
                                 correct <- rep(0,length(ordersInThere))
                               }
                               if (identical(order,NA)){
                                 order <- rep(NA, length(ordersInThere))
                               }
                               for (counter in 1:length(ordersInThere)){
                                 self$addTable(theData %>% filter(order == ordersInThere[counter]), dataMode = TRUE, whichInfo = whichInfo,
                                                correct = correct[counter], order = order[counter], rowNames = rowNames)
                               }
                             }
                           }
                         },
                         # finds and returns all places where the string toFind is present in the provided sequence
                         # order is important to combine with object table. Use type = ...  to give 'unique' identifier
                         # note that 'MOD_RES', 'CHAIN', 'COVER', etc are more or less reserved
                         # if !add, then returns the new table, if add, then the table gets added immeidately
                         # to the current (object) table
                         seqPart = function(theSequence = "", toFind =  "", type = "", description = "",
                                             order = 1, accession = "", entryName = "", taxid = 0, add = FALSE){
                           temp <- as.data.frame(str_locate_all(theSequence,toFind)[[1]]) %>%
                             arrange(start) %>% 
                             mutate(type = type,description = description, length = nchar(toFind), order = order,
                                    start = start, end = end + 1, accession = accession, entryName = entryName,
                                    taxid = taxid) %>%
                             rename(begin = start)
                           if (!add){
                             return(bind_rows(createTheoryTable(),temp))
                           } else {
                             self$addTable(bind_rows(createTheoryTable(),temp))
                           }
                         },
                         # essentially the function seqPart() applied to table (sequences have to be provided separately (1 per order!!)) 
                         seqParts = function(sequences = NA, toFind = "", orders = NA, type = "", description = "", add = FALSE){
                           if (identical(orders,NA)){
                             orders <- unique(self$table$order)
                           }
                           keepTable <- createTheoryTable()
                           for (counter in 1:(length(orders))){
                             tempTable <- (self$table %>% filter(order == orders[counter]))[1,]
                             keepTable <- bind_rows(keepTable,
                                                     self$seqPart(sequences[counter], toFind = toFind, type = type,
                                                                  description = description, order = orders[counter], 
                                                                  accession = tempTable$accession,
                                                                  entryName = tempTable$entryName,
                                                                  taxid = tempTable$taxid,
                                                                  add = FALSE
                                                                  )
                                           )
                           }
                           if (!add){
                             return(keepTable)
                           } else {
                             self$table <- bind_rows(self$table,keepTable)
                           }
                         },
                         # this wrapper uses the original draw_chains, but is prepared for experimentTables, which do not have entryNames
                         # the rows with type == "CHAINS" are expected to have the label to be used (eg name of the protein) in the
                         # description column
                         draw_chains = function(p, data = self$table, ...){
                           theChains <- which(is.na(data[data$type == "CHAIN",]$entryName))
                           if (length(theChains) > 0){
                             data[data$type == "CHAIN",]$entryName[theChains] <- data[data$type == "CHAIN",]$description[theChains]
                           }
                           return(drawProteins::draw_chains(p, data, ...))
                         }
                       ),
                       active = list(
                         # returns a list of all rows with type == "CHAIN"
                         proteins = function(value) {
                           if (missing(value)) {
                             return(self$table %>%
                                      filter(type == "CHAIN") %>%
                                      select(order, description, accession, entryName, taxid) %>%
                                      rename(name = description, entry = entryName, taxonomy = taxid))
                           } else {
                             # do nothing, read only
                           }
                         }
                         )
)

experimentChain <- R6Class("experimentChain",
                           inherit = chainData,
                           public = list(
                             createTable = function(){
                               self$table <- createExperimentTable()
                             },
                             # to add data to the experimentTable, care needs to be used with regard to the order column!
                             # data needs to be in the format produced by createExperimentTable():
                             # type = "COVER": to indicate the data (sequence, start aa, etc) relates to a part of protein identified
                             #                 a row with of this type contains info on modification (if present) in the mods column
                             #        "MOD_RES": to indicate the position and type of a modification
                             #        Note: other values are ignored
                             # description = in case type == "MOD_RES", then contains the name of the modification, begin == end, length == 1
                             #                       type == "COVER", then description is not used (usually contains name protein),
                             #                                        begin, end & length deal with the identified region of the protein
                             # begin = start or position of identified sequence or modification respectively
                             # end = end of identified sequence (otherwise not used)
                             # length = end of identified sequence (otherwise not used)
                             # order = see theoryTable, it is advisable to keep the order in all 3 tables the same!
                             # sequence = in case of type == "COVER" it will contain the amino acid sequence of the identified part of protein
                             # mods = in case of type == "COVER" it will contain the description of the modifications (if present) in the
                             #                                   identified part of protein. Ignored in all other cases
                             # expMr = experimental neutral mass of the ion
                             # calcMr = calculated (theoretical) neutral mass of the sequence (including modifications)
                             # rt = retention time of the identified ions
                             # expect = value relating to quality of identification (in case of Mascot data, this is actual expect value)
                             # varModPos = character string used to idnicate the position of modifications in a sequence, this comes from Mascot
                             #             originally and is non-essential for the functioning of the functions here
                             # a function that allows addition of multiple data.frames of data seems unnecessary at this moment
                             # arguments: order, if not NA, then the value here will be used to set/overwrite the order column in theData
                             addTable = function(theData = NA, order = NA){
                               if (!identical(theData,NA) & !identical(theData,NULL)){
                                 if (!identical(order,NA)){
                                   theData$order <- order
                                 }
                                 self$table <- bind_rows(self$table,theData)
                               } 
                             },
                             # same as addTheory, but now can work with either theData = list of filenames to load (dataMode = FALSE)
                             # or a data.frame containing more than one value for order
                             # note: order = list of numbers to use for the order of the different files/data
                             addTables = function(theData = NA, order = NA){
                               if (!identical(theData, NULL)){
                                 ordersInThere <- unique(theData$order)
                                 if (identical(order,NA)){
                                   order <- rep(NA, length(ordersInThere))
                                 }
                                 for (counter in 1:length(ordersInThere)){
                                   self$addTable(theData %>% filter(order == ordersInThere[counter]), order = order[counter])
                                 }
                               }
                             },
                             # finds and returns all places where the string toFind is present in the provided sequence
                             # order is important to combine with object table. Use type = ...  to give 'unique' identifier
                             # note that 'MOD_RES', 'CHAIN', 'COVER', etc are more or less reserved
                             # if !add, then returns the new table, if add, then the table gets added immeidately
                             # to the current (object) table
                             seqPart = function(theSequence = "", toFind =  "", type = "", description = "",
                                                 order = 1, accession = "", entryName = "", taxid = 0, add = FALSE){
                               temp <- as.data.frame(str_locate_all(theSequence,toFind)[[1]]) %>%
                                 arrange(start) %>% 
                                 mutate(type = type,description = description,length = nchar(toFind), order = order,
                                        sequence = toFind, start = start, end = end + 1) %>%
                                 rename(begin = start)
                               if (!add){
                                 return(bind_rows(createExperimentTable(),temp))
                               } else {
                                 self$addTable(bind_rows(createExperimentTable(),temp))
                               }
                             },
                             # essentially the function seqPart() applied to table (sequences have to be provided separately (1 per order!!)) or
                             # every order needs to have a chain+sequence row
                             seqParts = function(toFind = "", orders = NA, type = "", description = "", sequences = NA, add = FALSE){
                               if (identical(orders,NA)){
                                 orders <- unique(self$table$order)
                               }
                               tempTable <- createExperimentTable()
                               for (counter in 1:(length(orders))){
                                 tempTable  <- bind_rows(tempTable,
                                                         self$seqPart((self$table %>%
                                                                         filter(order == orders[counter], type == "CHAIN"))$sequence, # one chain per order!!
                                                                      toFind = toFind, type = type,
                                                                      description = description,
                                                                      order = orders[counter])
                                                         )
                               }
                               if (!add){
                                 return(tempTable)
                               } else {
                                 self$table <- bind_rows(self$table, tempTable)
                               }
                             },
                             # calculates coverage of protein sequence (which parts have been identified)
                             # dependent on: - a row with type = "CHAIN" & description = name protein & sequence = amino acid
                             #                 sequence protein (usually from database)
                             #               - column "type" containing entries "COVER" (usually) or as defined by the argument
                             #                 coverString
                             #               - colums "begin" and "end" which are used to determine identified parts of the sequence 
                             coverage = function(order = NA, coverString = "COVER"){
                               coverageTable <- data.frame(protein = as.character(), coverage = as.numeric(),
                                                           identifications = as.numeric(), stringsAsFactors = FALSE)
                               if (identical(order,NA)){
                                 ordersInThere <- unique(self$table$order)
                               } else {
                                 ordersInThere <- order
                               }
                               for (counter0 in 1:(length(ordersInThere))){
                                 lengthProt <- self$table %>%
                                   filter((type == "CHAIN") & (order == ordersInThere[counter0]))
                                 if (nrow(lengthProt)>0){                                # else : current order has no chain
                                   coverage <- rep(FALSE,nchar(lengthProt$sequence[1]))  # [1] in case there's more than one chain present
                                   tempTable <-  self$table %>% 
                                     filter((order == ordersInThere[counter0]) & (type == coverString)) %>%
                                     select(begin, end)
                                   if (nrow(tempTable)>0){
                                     for (counter in 1:nrow(tempTable)){
                                       for (counter2 in tempTable[counter,]$begin:(tempTable[counter,]$end-1)){
                                         coverage[counter2] <- TRUE
                                       }
                                     }
                                     coverageTable <- bind_rows(coverageTable,
                                                                data.frame(protein = (self$table %>%
                                                                                        filter((type == "CHAIN") &
                                                                                                 (order == ordersInThere[counter0])))$description,
                                                                           coverage = (sum(coverage)/length(coverage))*100,
                                                                           identifications = nrow(tempTable))
                                     )
                                   } else {
                                     coverageTable <- bind_rows(coverageTable,
                                                                data.frame(protein = (self$table %>%
                                                                                        filter((type == "CHAIN") &
                                                                                                 (order == ordersInThere[counter0])))$description,
                                                                           coverage = NA,
                                                                           identifications = NA)
                                     )
                                   }
                                 }
                               }
                             return(coverageTable)
                             },
                             # this wrapper uses the original draw_chains
                             # the rows with type == "CHAINS" are expected to have the label to be used (eg name of the protein) in the
                             # description column
                             draw_chains = function(p, data = self$table, ...){
                               data$entryName <- ""
                               data[data$type == "CHAIN",]$entryName <- data[data$type == "CHAIN",]$description
                               return(drawProteins::draw_chains(p, data, ...))
                             }
                           ),
                           active = list(
                             # returns a list of all rows with type == "CHAIN"
                             proteins = function(value) {
                               if (missing(value)) {
                                 return(self$table %>%
                                          filter(type == "CHAIN") %>%
                                          select(order, description, sequence) %>%
                                          rename(name = description))
                               } else {
                                 # do nothing, read only
                               }
                             }
                           )
)

# a function that takes a series of colors/fills/shapes and adds it as a legend to the ggplot object p
# note: each row of the legend data.frame is an item in the legend to be created, so the color, shape, etc
#       of a single row 'belong' together
# note: the row-order determines the order in the legend
# note: use of this function will result in a series of warning messages (one for every item in the legend)
#       this is 'normal' because the function uses a trick to do its thing. Use suprressWarnings(print())
#       to prevent seeing the messages while running code. In R markdown use the option 'warning=FALSE'
#       to not see the messages
customLegend <- function(p, legend = data.frame(labels = NA, colors = NA, fills = NA, shapes = NA, sizes = NA),
                         legend.title = "Legend", legend.title.face = "bold", legend.title.size = 12,
                         legend.element.size = 12, legend.position = "bottom"){
  for (counter in 1:(nrow(legend))){
    fakedf <- data.frame(x1x = 0,
                         y1y = 0,
                         group = legend$labels[counter])
    p <- p + geom_point(data = fakedf, aes(x = x1x, y = y1y,
                                           color = group,
                                           fill = group,
                                           shape = group,
                                           size = group))
  }
  p <- p +
    scale_color_manual(values = legend$colors, breaks = legend$labels) +
    scale_fill_manual(values = legend$fills, breaks = legend$labels) +
    scale_shape_manual(values = legend$shapes, breaks = legend$labels) +
    scale_size_manual(values = legend$sizes, breaks = legend$labels) + 
    labs(color = legend.title, fill = legend.title, shape = legend.title, size = legend.title) +
    theme(legend.title=element_text(size = legend.title.size, face = legend.title.face),
          legend.margin=margin(l=0),
          legend.text=element_text(size = legend.element.size),
          legend.position = legend.position)
  return(p)
}



