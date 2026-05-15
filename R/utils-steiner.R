#' Function to aproximate a Steiner tree based on the KB heuristic approach
#'
#' See: Afshin Sadeghi and Holger Froehlich, "Steiner tree
#'  methods for optimal sub-network identification: an empirical study",
#'  BMC Bioinformatics 2013 14:144 doi:10.1186/1471-2105-14-144.
#'
#' @noRd
steinertree_KB <- function(terminals, g) {

  # Helper function:
  # Given a graph and a set of nodes, return a list of connected components
  get_connected_comps <- function(g, all_nodes) {
    g_induced <- induced_subgraph(g, vids = all_nodes)
    comps <- components(g_induced)
    lapply(1:comps$no, function(i, comps) {names(comps$membership)[comps$membership == i]}, comps)
  }

  # Make a Steiner Tree from every terminal
  subtrees <- get_connected_comps(g, terminals)

  # Will hold all terminals not in a steiner subtree
  # (Note: This is a "trick" for shortening time on searching for shortest paths.
  #  It is equivalent to the original problem)
  nsubtrees <- lapply(subtrees, function (r) setdiff(terminals, r))

  # While no subtree includes all terminals...
  while (length(subtrees) > 1) {

    # Find shortest paths between different Steiner Trees
    paths <- lapply(
      1:length(subtrees),
      function (r) lapply(
        intersect(subtrees[[r]],terminals),
        function (x, y, g) { all_shortest_paths(g, from = x, to = y)$res },
        y = nsubtrees[[r]],
        g = g
      ) %>% unlist(recursive=FALSE)
    )

    # Compute the length of each shortest path
    t <- sapply(1:length(paths), function (r) sapply(paths[[r]], length))

    # Compute a minimum for each set of lengths from each Steiner tree to other trees
    t2 <- sapply(t, function (x) min(x))
    shortest_paths_len <- min(t2)

    # Iterate over paths. If path length == shortest_paths_len, return the corresponding nodes.
    nodes_on_shortest_paths <- lapply(
      paths,
      function (paths_i, shortest_paths_len) {
        Filter(
          Negate(is.null),
          lapply(
            paths_i,
            function(paths_i_j, shortest_paths_len) {
              if (length(paths_i_j) == shortest_paths_len)
                return(paths_i_j$name)
              return(c())
            },
            shortest_paths_len
          ) # %>% unlist()
        )
      },
      shortest_paths_len
    ) %>% unlist(recursive = F)

    # Randomly (!) select only one of all paths identified to add to the tree
    np <- length(nodes_on_shortest_paths)
    nodes_on_shortest_paths <- nodes_on_shortest_paths[[sample(np, size = 1)]]

    # plot(induced_subgraph(g, unique(c(nodes_on_shortest_paths, unlist(subtrees)))))
    subtrees <- get_connected_comps(g, unique(c(nodes_on_shortest_paths, unlist(subtrees))))
    nsubtrees <- lapply(subtrees, function (r) setdiff(terminals, r))
  }

  # Optimize
  stree <- mst(induced_subgraph(g, subtrees[[1]]))
  non_temrinal_leafs <- igraph::degree(stree, v = setdiff(V(stree)$name, terminals))
  non_temrinal_leafs <- names(non_temrinal_leafs)[non_temrinal_leafs == 1]
  while(length(non_temrinal_leafs) > 0) {
    stree <- delete_vertices(stree, v = non_temrinal_leafs)
    non_temrinal_leafs <- igraph::degree(stree, v = setdiff(V(stree)$name, terminals))
    non_temrinal_leafs <- names(non_temrinal_leafs)[non_temrinal_leafs == 1]
  }

  return(stree)
}

#' Recieves a list of modules, each including a set of nodes, and computes a
#'   steiner tree to assure each module forms a connected subgraph
#'
#' @noRd
complete_modules_with_steiner <- function(
    g,
    modules_overview,
    module_assignments
    ) {
  # Place holder
  complete_modules <- data.frame()

  # Iterate over modules
  for (mod_id in modules_overview$module_id) {
    cat('.')
    if (mod_id%%10 == 0) cat('\n')
    
    # Fetch module seeds (red nodes)
    module_anchor_nodes <- module_assignments %>%
      filter(module_id == mod_id) %>%
      pull(name) %>%
      as.character()

    # Use them as terminal nodes for the steiner tree algorithm
    stree <- steinertree_KB(
      terminals = module_anchor_nodes,
      g = g
    )
    module_all_nodes <- V(stree)$name

    # Save full modules
    complete_modules <- bind_rows(
      complete_modules,
      data.frame(
        node = module_all_nodes,
        is_anchor = V(g)$anchor[match(module_all_nodes, V(g)$name)] == 1,
        module_id = mod_id
      )
    )
  }

  # # Some strict clustering methods seem to 'cut' larger modules, which after 
  # #. completion with Steiner trees get stitched back together. 
  # #. We identify these stitched modules and update their module IDs. 
  # 
  # ## First identify pairs of modules sharing at least one node
  # module_edges <- complete_modules %>%
  #   select(node, module_id) %>%
  #   inner_join(complete_modules, by = "node", relationship = "many-to-many") %>%
  #   filter(module_id.x < module_id.y) %>%
  #   distinct(module_id.x, module_id.y)
  # 
  # ## Build a graph of module overlaps and extract connected components (= modules to merge)
  # g_tmp <- graph_from_data_frame(
  #   module_edges,
  #   directed = FALSE,
  #   vertices = tibble(module_id = unique(complete_modules$module_id))
  # )
  # module_id_old_to_new <- tibble(
  #   module_id = as.integer(names(components(g_tmp)$membership)),
  #   new_module_id = components(g_tmp)$membership
  # )
  # ## Attach new module IDs back to original table
  # complete_modules <- complete_modules %>%
  #   left_join(module_id_old_to_new, by = "module_id")
  
  # Add p-values and other info
  complete_modules$pval <- vertex_attr(g, 'pval', complete_modules$node)
  complete_modules$type <- vertex_attr(g, 'type', complete_modules$node)

  return(complete_modules)
}
