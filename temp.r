library("Rook")


app <- function(env) {

  req <- Request$new(env)
  res <- Response$new()

  input = fromJSON(rawToChar(env[['rook.input']]$postBody))

  upgenes <- input$upgenes
  downgenes <- input$downgenes
  direction <- input$direction
  species <- input$species
  
  gene <- input$gene
  sigtype <- input$type 
  
  signatureInput <- input$signature
  siggenes <- input$siggenes
  
  if(!is.null(gene)){
    print("gene correlation")
    ll = list()
    
    if(gene %in% rownames(cc)){
      ll[[1]] = rownames(cc)
      ll[[2]] = cc[gene,]
    } else if(gene %in% rownames(human_correlation)){
      ll[[1]] = rownames(human_correlation)
      ll[[2]] = human_correlation[gene,]
    } else{
      ll[[1]] = ""
    }
    
    body <- toJSON(ll)
    ret <- paste0(body, "\n");
    
    res$header("Content-type", "application/json")
    res$header("Access-Control-Allow-Origin", "*")
    res$write(ret)
    res$finish()
    
  } else if(sigtype == "full_signature"){
    
    if(species == "human"){
      genes = colnames(human_transform)
      expression = jl_human_expression
      transform = human_transform
      rmean = rmean_human
      rsd = rsd_human
      reference_dist = reference_dist_human
    } else if(species == "mouse"){
      genes = colnames(mouse_transform)
      expression = jl_mouse_expression
      transform = mouse_transform
      rmean = rmean_mouse
      rsd = rsd_mouse
      reference_dist = reference_dist_mouse
    }
    
    signature = rep(0, length(genes))
    names(signature) = genes
    
    intergene = which(siggenes %in% genes)
    signature[siggenes[intergene]] = as.numeric(signatureInput[intergene])
    
    r2 = rank(signature)
    vec = (reference_dist[r2] - rmean)/rsd
    
    jl_vec = transform %*% vec
    
    similarity = scale(c(cor(jl_vec, expression)))
    names(similarity) = colnames(expression)
    
    ll = list()
    ll[[1]] = input$signatureName
    names(ll)[1] = "name"
    
    ll[[2]] = as.numeric(gsub("GSM", "", names(similarity)))
    names(ll)[2] = "samples"
    
    ll[[3]] = similarity
    names(ll)[3] = "similarity"
    
    body <- toJSON(ll)
    ret <- paste0(body, "\n");
    
    res$header("Content-type", "application/json")
    res$write(ret);
    
    res$finish()
    
  } else if(sigtype == "geneset"){
    if(species == "human"){
      genes = colnames(human_transform)
      expression = jl_human_expression
      transform = human_transform
    } else if(species == "mouse"){
      genes = colnames(mouse_transform)
      expression = jl_mouse_expression
      transform = mouse_transform
    }
    
    vec = genes %in% upgenes
    vec[genes %in% downgenes] = -1
    names(vec) = genes
    jl_vec = transform %*% vec
    
    similarity = scale(c(cor(jl_vec, expression)))
    names(similarity) = colnames(expression)
    ww = which(similarity > 2.5)
    
    samples = names(similarity)[ww]
    
    ll = list()
    
    ll[[1]] = input$signatureName
    names(ll)[1] = "name"
    
    ll[[2]] = as.numeric(gsub("GSM","",samples))
    names(ll)[2] = "samples"
    
    body <- toJSON(ll)
    ret <- paste0(body, "\n");
    
    res$header("Content-type", "application/json")
    
    res$write(ret);

    res$finish()
  }
}



