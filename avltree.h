#include<iostream>
using namespace std;
 
 
typedef size_t ElenmentType; 
typedef struct AVLNode{
	int depth;
    int tombstone; // 0:delete 1:insert
	struct AVLNode *left;
	struct AVLNode *right;
	struct AVLNode *parent;
	ElenmentType value; 
}AVLtree,Tree;

Tree* avlInit(ElenmentType value){
	Tree* root=new Tree();
	root->parent = NULL;
    root->depth = 0;
    root->left = root->right = NULL;
    root->value=value; 
    root->tombstone = 1; // insertion
    return root;
}

Tree* avlInit(ElenmentType value, int tombstone){
	Tree* root=new Tree();
	root->parent = NULL;
    root->depth = 0;
    root->left = root->right = NULL;
    root->value=value; 
    root->tombstone = tombstone; 
    return root;
}
 
Tree* LL_rotate(Tree *root){
    Tree *temp;
	temp = root->left;
	root->left = temp->right;
	temp->right = root;
	return temp; 
   
}

Tree* RR_rotate(Tree * root){
    Tree* temp;
    temp = root->right;
    root->right = temp->left;
    temp->left = root;
    return temp;
}

Tree* LR_rotate(Tree* root){
	Tree* temp;
	temp = root->left;
    root->left =RR_rotate(temp);
    return LL_rotate(root);
} 

Tree* RL_rotate(Tree* root){
	 Tree* temp;
	 temp = root->right;
    root->right=LL_rotate(temp);
    return RR_rotate(root);
} 
 
int height(Tree* root)
{
    if (root == NULL)
        return 0;
    int max;
    int left=height(root->left);
    int right=height(root->right);
    if(left>=right)
    max=left;
    else
    max=right;
    return max+1;
}

int  GetSumOfLeafNode(Tree* root)
{
	if(root == NULL)
		return 0;
	
	if(root->left == NULL && root->right == NULL)
		return 1;
	else
	{
		return GetSumOfLeafNode(root->left) 
			+ GetSumOfLeafNode(root->right);
	}
}

int diff(Tree* root)
{
    return height(root->left) - height(root->right);
}
 
Tree* avlBalance(Tree* root)
{
    int balanceFactor = diff(root);
    if (balanceFactor > 1)
    {
        if (diff(root->left) > 0)
		
            root=LL_rotate(root);
        else
		
            root=LR_rotate(root);
    }
    else if (balanceFactor < -1)
    {
        if (diff(root->right) > 0)
		
            root=RL_rotate(root);
        else
		
            root=RR_rotate(root);
    }
    return root;
} 

Tree* avlInsert(Tree* root, ElenmentType k, int tombstone)
{
    if (NULL == root)
    {
        root = avlInit(k, tombstone);
        if(root==NULL)
        	cout<<"creat failed"<<endl; 
        return root;
    }
    else if (k < root->value)
    {
        root->left = avlInsert(root->left, k, tombstone);
        root = avlBalance(root);
    }
    else if (k>root->value)
    {
        root->right = avlInsert(root->right, k, tombstone);
        root = avlBalance(root);
    }
    return root;
} 

void displayTree_Pre(Tree* node){
        if(node == NULL) return;
        cout<<node->value<<" ";
        if(node->left != NULL){
            displayTree_Pre(node->left);
        }
        if(node->right != NULL){
            displayTree_Pre(node->right);
        }
}

void displayTree_Mid(Tree* node){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Mid(node->left);
        }
        cout<<node->value<<" ";
        if(node->right != NULL){
            displayTree_Mid(node->right);
        }
}


void displayTree_Mid(Tree* node, vector<size_t> &vec, size_t src){
    if(node == NULL) return;    
    if(node->left != NULL){
        displayTree_Mid(node->left, vec, src);
    }

    // if(src == 9858)
    //     cout<<node->value<<" ";

    if (vec.size() >= 1) 
    {
        if ( (vec[vec.size()-1] != src) && (node->value - vec[vec.size()-1] > 1) )
        {
            vec.push_back(src);  
        }
    }
    vec.push_back(node->value);
    if(node->right != NULL){
        displayTree_Mid(node->right, vec, src);
    }
}


void displayTree_Mid(Tree* node, size_t src, size_t &edge_num){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Mid(node->left, src, edge_num);
        }

        edge_num++;

        if(node->right != NULL){
            displayTree_Mid(node->right, src, edge_num);
        }
}

void displayTree_Post(Tree* node){
        if(node == NULL) return;
        if(node->left != NULL){
            displayTree_Post(node->left);
        }
        if(node->right != NULL){
            displayTree_Post(node->right);
        }
        cout<<node->value<<" ";
}

Tree* binaryTreeSearch(Tree *node,int value){
	if(node->value==value)
		return node;
	
	else if(node->value>value){
		if(node->left!=NULL)
			return binaryTreeSearch(node->left,value);
		else return NULL;
	}
	
	else{
		if(node->right!=NULL)
			return binaryTreeSearch(node->right,value);
		else
			return NULL;
	}
}
 

ElenmentType tree_max(Tree *node){
	int value;
	value=node->value;
	if(node->right!=NULL)
	return tree_max(node->right);
	else
	return value;
	
}
 
ElenmentType tree_min(Tree *node){
	int value;
	value=node->value;
	if(node->left!=NULL)
	return tree_min(node->left);
	else
	return value;
}
 

Tree* avlDelete(Tree *root, const ElenmentType k)
{
    if (NULL == root)
        return root;
    if (!binaryTreeSearch(root,k))
    {
        cout<<"Delete failed"<<endl;
        return root;
    }
 
    if (k == root->value)
    {
        if (root->left!=NULL&&root->right!=NULL)
        {
            if (diff(root) > 0)
            {
                root->value = tree_max(root->left);
                root->left = avlDelete(root->left, root->value);
            }
            else
            {
                root->value = tree_min(root->right);
                root->right = avlDelete(root->right, root->value);
            }
        }
        else
        {
                Tree * tmp = root;
                root = (root->left) ? (root->left) :( root->right);
                delete tmp;
                tmp = NULL;
        }
    }
    
    else if (k < root->value)
    {
        root->left = avlDelete(root->left, k);
        
		if (diff(root) < -1)
        {
            if (diff(root->right) > 0)
            {
                root = RL_rotate(root);
            }
            else
            {
                root = RR_rotate(root);
            }
        }
    }
    
    else
    {
        root->right = avlDelete(root->right, k);
        
        if (diff(root) > 1)
        {
            if (diff(root->left) < 0)
            {
                root = LR_rotate(root);
            }
            else
            {
                root = LL_rotate(root);
            }
        }
    }
    return root;
}
 
 
 
// int main(){
// 	int a[10]={3,5,4,2,6,8,10,1,7,9};
// 	Tree *root=NULL;
// 	for(int i=0;i<10;i++){
// 		root = Insert(root,a[i]);
// 	}
// 	cout<<"Preorder Traversal: ";
// 	displayTree_Pre(root);
// 	cout<<endl; 
// 	cout<<"Midorder Traversal: ";
// 	displayTree_Mid(root);
// 	cout<<endl; 
// 	cout<<"Tree Height:"<<height(root)<<endl; 
// 	int value = 33;
// 	cout<<"Max key: "<<tree_max(root)<<endl;
// 	cout<<"Min key: "<<tree_min(root)<<endl;
// 	cout<<"Num of leaf: "<<GetSumOfLeafNode(root)<<endl; 
// 	Tree* obj;
// 	if((obj=binaryTreeSearch(root,value))==NULL){
// 		cout<<value<<"No"<<endl;
// 	}
// 	else 
// 	    cout<<value<<"Yes"<<endl;
// 	root = Insert(root,5);
// 	cout<<"insert 5, Preorder Traversal: ";
//     displayTree_Pre(root);
// 	cout<<endl;
// 	cout<<"insert 5, Midorder Traversal: ";
// 	displayTree_Mid(root);
// 	cout<<endl;
// 	int del=79;
// 	if (!binaryTreeSearch(root,del))
//     {
//         cout<<"Delete failed, no result"<<endl;
//     }else{
//     	Delete(root,del);
//     	cout<<"Delete node, Preorder Traversal: ";
//     	displayTree_Pre(root);
//     	cout<<endl;
//     	cout<<"Delete node, Midorder Traversal: ";
//     	displayTree_Mid(root);
//     	cout<<endl;
// 	} 
// 	return 0;
// }
