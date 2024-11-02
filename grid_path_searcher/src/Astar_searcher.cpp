#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    //lower
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    //upper
    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    //创建每个网格
    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    //创建容器储存每个网格的节点状态
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                //离散网格向物理坐标转换
                Vector3d pos = gridIndex2coord(tmpIdx);
                //储存node的网格index值和实际坐标值
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLX_SIZE ; i++)
        for(int j=0; j < GLY_SIZE ; j++)
            for(int k=0; k < GLZ_SIZE ; k++)
                resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                //if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    //+0.5表示取grid的中心处
    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();
    /*
    *
    STEP 4: finish AstarPathFinder::AstarGetSucc yourself 
    please write your code below
    *
    *
    */

    //explore 26 node in 3D scene
    Vector3i current_index = currentPtr->index;
    Vector3i index_lower, index_upper;

    index_lower << min(max(current_index(0)-1,0),GLX_SIZE-1), min(max(current_index(1)-1,0),GLY_SIZE-1), min(max(current_index(2)-1,0),GLZ_SIZE-1);
    index_upper << min(max(current_index(0)+1,0),GLX_SIZE-1), min(max(current_index(1)+1,0),GLY_SIZE-1), min(max(current_index(2)+1,0),GLZ_SIZE-1);

    // index_lower << current_index(0)-1, current_index(1)-1, current_index(2)-1;
    // index_upper << current_index(0)+1, current_index(1)+1, current_index(2)+1;

    for(int i = index_lower(0); i <= index_upper(0); i++){
        for(int j = index_lower(1); j <= index_upper(1); j++){
            for(int k = index_lower(2); k <= index_upper(2); k++)
            {
                Vector3i neighbor_index(i,j,k);
                if(isFree(neighbor_index) && neighbor_index != current_index){
                    neighborPtrSets.emplace_back(GridNodeMap[i][j][k]);
                    //Tie breaker
                    //edgeCostSets.emplace_back(sqrt(pow(i-current_index(0),2)+pow(j-current_index(1),2)+pow(k-current_index(2),2))*1.0001);
                    edgeCostSets.emplace_back(sqrt(pow(i-current_index(0),2)+pow(j-current_index(1),2)+pow(k-current_index(2),2)));
                }
            }
        }
    }

}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    /* 
    choose possible heuristic function you want
    Manhattan, Euclidean, Diagonal, or 0 (Dijkstra)
    Remember tie_breaker learned in lecture, add it here ?
    *
    *
    *
    STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    please write your code below
    *
    *
    */
    //欧式距离
    // double heuristic_distance_ = (node1->coord - node2->coord).norm();
    // return heuristic_distance_;

    //Octile Distance
    Eigen::Vector3d node1_coordinate = node1->coord;
    Eigen::Vector3d node2_coordinate = node2->coord;

    double dx = fabs(node1_coordinate[0] - node2_coordinate[0]);
    double dy = fabs(node1_coordinate[1] - node2_coordinate[1]);
    double dz = fabs(node1_coordinate[2] - node2_coordinate[2]);

    double D1 = 1;
    double D2 = sqrt(2);
    double D3 = sqrt(3);

    double Diagonal_min = min(min(dx, dy), dz);
    double Diagonal_max = max(max(dx, dy), dz);
    double Diagonla_mid = dx + dy + dz - Diagonal_max - Diagonal_min;

    double epsilon = 1.001;  // slightly favor nodes closer to the goal
    double heuristic_distance_ = ((D3 - D2) * Diagonal_min + (D2 - D1) * Diagonla_mid + Diagonal_max) * epsilon;

    //double heuristic_distance_ = (D3-D2)*Diagonal_min + (D2-D1)*Diagonla_mid + Diagonal_max;

    //ROS_INFO("heu:%lf",heuristic_distance_);
    return heuristic_distance_;

}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    

    //index of start_point and end_point
    //将实际距离转化为网格地图
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;

    //position of start_point and end_point
    //实际地图转化为栅格地图
    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    //创建新gridnode
    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    // currentPtr represents the node with lowest f(n) in the open_list
    //创建当前节点与当前节点的邻居节点
    GridNodePtr currentPtr  = NULL;
    GridNodePtr neighborPtr = NULL;

    //初始化节点内的信息，包括路径代价与启发信息以及当前节点在openlist或closelist
    //put start node in open set
    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr,endPtr);   
    //STEP 1: finish the AstarPathFinder::getHeu , which is the heuristic function
    startPtr -> id = 1; 
    startPtr -> coord = start_pt;
    //创建键值对储存在openlist中，键值对有fscore和节点的指针组成
    openSet.insert( make_pair(startPtr -> fScore, startPtr) );
    /*
    *
    STEP 2 :  some else preparatory works which should be done before while loop
    please write your code below
    *
    *
    */
    //邻居节点集合与代价集合
    vector<GridNodePtr> neighborPtrSets;
    vector<double> edgeCostSets;

    // this is the main loop
    while ( !openSet.empty() ){
        /*
        *
        *
        step 3: Remove the node with lowest cost function from open set to closed set
        please write your code below
        
        IMPORTANT NOTE!!!
        This part you should use the C++ STL: multimap, more details can be find in Homework description
        *
        *
        */
        int time;
        time ++;
        ROS_INFO("TIME:%d",time);
        //由于openset中是fscore和节点指针的键值对，begin() 返回 openSet 的第一个元素的迭代器，使用begin自动排序出fscore最小的节点指针，取第二个参数即取出节点指针
        currentPtr = openSet.begin()->second;
        //openSet.begin()返回的是代价最低的那个节点,即将该节点从openlist中移出
        openSet.erase(openSet.begin());
        //add to close set
        currentPtr->id = -1;
        ROS_INFO("MINFS:%f",currentPtr->fScore);

        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
            return;
        }
        //get the succetion 寻找更新邻居节点
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);  //STEP 4: finish AstarPathFinder::AstarGetSucc yourself     

        /*
        *
        *
        STEP 5:  For all unexpanded neigbors "m" of node "n", please finish this for loop
        please write your code below
        *        
        */         
        //在所有邻居节点中筛选，符合的加入openlist，同时更新fscore值
        for(int i = 0; i < (int)neighborPtrSets.size(); i++){

            neighborPtr = neighborPtrSets[i];
            /*
            *
            *
            Judge if the neigbors have been expanded
            please write your code below
            
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : unexpanded
            neighborPtrSets[i]->id = 1 : expanded, equal to this node is in close set
            *        
            */
            //第一次探索的节点，之前未发现过，更新参数并加入openlist
            if(neighborPtr -> id == 0){ //discover a new node, which is not in the closed set and open set
                /*
                *
                *
                STEP 6:  As for a new node, do what you need do ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
                neighborPtr->cameFrom = currentPtr;
                neighborPtr->gScore = currentPtr->gScore + edgeCostSets[i];
                neighborPtr->fScore = neighborPtr->gScore + getHeu(neighborPtr,endPtr);
                neighborPtr->id = 1;
                openSet.insert( make_pair(neighborPtr -> fScore, neighborPtr) );
                //continue;
            }
            //之前探索过的节点，若经过当前路径的总代价小于之前探索结果的代价，则更新代价值，否则跳过
            else if(neighborPtr -> id == 1){ //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding
                /*
                *
                *
                STEP 7:  As for a node in open set, update it , maintain the openset ,and then put neighbor in open set and record it
                please write your code below
                *        
                */
                if(neighborPtr->gScore > currentPtr->gScore + edgeCostSets[i])
                {
                    neighborPtr->cameFrom = currentPtr;
                    neighborPtr->gScore = currentPtr->gScore + edgeCostSets[i];
                    neighborPtr->fScore = neighborPtr->gScore + getHeu(neighborPtr,endPtr);
                }
                //continue;
            }
            else{//this node is in closed set
                /*
                *
                please write your code below
                *        
                */
                continue;
            }
            ROS_INFO("FSCORE:%f",neighborPtr->fScore);
        }      
    }
    //if search fails
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
    //储存路径点的物理坐标
    vector<Vector3d> path;
    //储存路径的节点
    vector<GridNodePtr> gridPath;
    /*
    *
    *
    STEP 8:  trace back from the curretnt nodePtr to get all nodes along the path
    please write your code below
    *      
    */

    //定位当前目标节点指针
    auto pptr = terminatePtr;
    //循环追溯父节点以找到完整路径的所有节点
    while(pptr->cameFrom != NULL){
        gridPath.push_back(pptr);
        pptr = pptr->cameFrom;
    }

    //提取路径节点中的坐标点
    for (auto ptr: gridPath)
        path.push_back(ptr->coord);
        
    //由于从终点往起点倒推的路径是反的，将坐标序列进行倒序排列得到从起点到终点的路径的坐标序列
    reverse(path.begin(),path.end());

    return path;
}